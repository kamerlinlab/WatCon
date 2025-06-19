from pymol import cmd

def save_secondary_structure(outfile="ss_output.txt", selection="name CA"):
    with open(outfile, "w") as f:
        def writer(resn, resi, ss):
            f.write(f"{resn}{resi:>4}  {ss}\n")
        cmd.iterate(selection, "writer(resn, resi, ss)", space={"writer": writer})

if __name__ == '__main__':
    save_secondary_structure()