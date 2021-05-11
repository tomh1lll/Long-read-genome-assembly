import sys

threads=sys.argv[1]
genome=sys.argv[2]
outputDir=sys.argv[3]
inputFiles=sys.argv[4]
out=open(sys.argv[5],"w")

out.write("[General]\n")
out.write("job_type = local\n")
out.write("job_prefix = nextPolish\n")
out.write("task = best\n")
out.write("rewrite = yes\n")
out.write("rerun = 3\n")
out.write("parallel_jobs = " + threads + "\n")
out.write("multithread_jobs = 4\n")
out.write("genome = " + genome + "\n")
out.write("genome_size = auto\n")
out.write("workdir = " + outputDir + "\n")
out.write("polish_options = -p {multithread_jobs}\n")
out.write("\n")
out.write("[lgs_option]\n")
out.write("lgs_fofn = " + inputFiles + "\n")
out.write("lgs_options = -min_read_len 1k -max_depth 200\n")
out.write("lgs_minimap2_options = -x map-pb\n")
