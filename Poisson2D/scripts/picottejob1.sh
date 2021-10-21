#!/bin/bash
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=grs53@drexel.edu
#SBATCH --account=math540Prj
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=48
#SBATCH --time=1:00:00
#SBATCH --partition=def

module load julia

export WORKDIR=/ifs/groups/math540Grp/simpson

mkdir -p ${WORKDIR}

printf 'Starting...\n'

julia -t 48 --project=../ solve1.jl

printf 'Moving data...\n'

mv *.jld2 ${WORKDIR}
