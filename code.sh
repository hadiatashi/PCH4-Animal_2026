#!/bin/bash
export LC_ALL=C
ulimit -s unlimited
export OMP_NUM_THREADS=12
export KMP_STACKSIZE=300G

# STEP 1: inb.inb (same for all traits, run once)
awk '
NR == FNR { inb[$1] = $2; next }
{ if ($9 in inb) print $1, inb[$9], $9 }
' ../renf90.inb ../pedi.1-3.ren > inb.inb

for t in $(seq 8 40); do
    tt=$(printf "%02d" "$t")
    echo "=== Trait $t ==="

    # Read this trait's coefficients
    read a b < <(awk -v t="$t" '$1 == t {print $2, $3}' coef.txt)
    if [ -z "$a" ] || [ -z "$b" ]; then
        echo "No coefficients for trait $t in coef.txt, skipping"; continue
    fi

    for f in ${t}_priors ${t}_varG ${t}_varR; do
        [ -f "$f" ] || { echo "Missing $f, skipping trait $t"; continue 2; }
    done

    # STEP 2: solutions
    awk '
    NR == 1 {
        for (i=1; i<=NF; i++) {
            if ($i == "trait") trait=i
            if ($i == "effect") effect=i
            if ($i == "level") level=i
            if ($i == "sol") sol=i
            if ($i == "REL") rel=i
            if ($i == "REL_MACE") rel_mace=i
        }
        print "trait effect level sol new_rel"
        next
    }
    {
        if ($rel_mace != "" && $rel_mace != 0) new_rel = $rel_mace
        else new_rel = $rel
        printf "%d %d %d %s %.5f\n", $trait, $effect, $level, $sol, new_rel
    }
    ' ${t}_priors > solutions

    # STEP 3: genomic REL
    /home/ULG/f031707/progs/genomicrel --ntrait 1 --effect 1 \
        --rel solutions --ped pedi.1-3.ren --inb inb.inb \
        --gvcv ${t}_varG --rvcv ${t}_varR grel \
        --gi Gi --xref code.snp_subset_XrefID

    if [ ! -f geno_rel.dat ]; then
        echo "geno_rel.dat not created for trait $t, skipping"; continue
    fi

    awk -v a="$a" -v b="$b" '
    NR == FNR {
        key = $1 "_" $2
        genrel[key] = a + b * $3
        next
    }
    FNR == 1 { print "trait effect level sol new_rel final_rel"; next }
    {
        key = $1 "_" $3
        new_rel = $5
        if (key in genrel) genomic_rel = genrel[key]
        else genomic_rel = 0
        final_rel = (new_rel > genomic_rel) ? new_rel : genomic_rel
        printf "%-2d %-2d %-8d %8.2f %6.2f %6.2f\n", $1, $2, $3, $4, new_rel, final_rel
    }
    ' geno_rel.dat solutions > updated_solutions_${tt}

    rm -f erc_add.dat erc_f.dat erc_geno.dat erc_geno_ind_1.dat erc_gg.dat \
          erc_ss.dat geno_rel.dat log_genomicrel.dat pev_geno_ind_1.dat \
          rel_g.dat solutions
done

echo "Finished."