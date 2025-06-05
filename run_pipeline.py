import os,sys,glob,argparse,subprocess

import argparse
def main():
    parser = argparse.ArgumentParser(description="Run full KEGG KO-based pathway and reaction inference pipeline.")
    parser.add_argument("--ko_dir", required=True, help="Directory containing .ko files for each genome.")
    parser.add_argument("--module_file", required=True, help="KEGG module definition file (e.g., ko_pathway.list).")
    parser.add_argument("--keyenzyme_dir", required=True, help="Directory containing KOs for simple reactions")
    parser.add_argument("--output_dir", required=True, help="Output directory for results.")
    args  = parser.parse_args()  

    if not os.path.exists(args.module_file): raise ValueError(f"No KEGG Database file for module definition was found from given path (args.module_file)")

    if not os.path.exists(args.output_dir): os.mkdir(args.output_dir)

    subprocess.run(['python', './scripts/1_track_module_steps.py', f"--ko_dir={args.ko_dir}", f"--module_file={args.module_file}", f"--output={args.output_dir}/Steps.txt"])
    subprocess.run(['python', './scripts/2_evalulate_module_completion.py', f"--input={args.output_dir}/Steps.txt", f"--output={args.output_dir}/Pathways.txt"])
    subprocess.run(['python', './scripts/3_infer_simple_reactions.py', f"--ko_dir={args.ko_dir}", f"--keyenzyme_dir={args.keyenzyme_dir}", f"--output={args.output_dir}/Reactions.txt"])
    

if __name__ == "__main__": main()

