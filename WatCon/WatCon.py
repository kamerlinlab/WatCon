import os, sys
import argparse

def parse_inputs(filename):
    kwargs = {}
    structure_type = None
    with open(filename, 'r') as FILE:
        for line in FILE:
            if line.startswith(';') or ':' not in line:
                continue
            elif 'structure_type' in line:
                structure_type = line.split(':')[1].split()[0]
            else:
                kw = line.split(':')[0]
                kw_value = line.split(':')[1].split()[0]
                kwargs[kw] = kw_value
    return structure_type, kwargs

def check_conditions(kwargs):
    print('Checking conditions from input file')
    #Consider making some text to check conditions
    pass
    
def run_watcon(structure_type, kwargs):
    kwargs['return_network'] = False #If input files are being used, don't return networks
    if structure_type == 'static':
        from generate_static_networks import initialize_network
    else:
        from generate_dynamic_networks import initialize_network
    results = initialize_network(*kwargs)
    print(results)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform analysis using WatCon')
    parser.add_argument('--input', type=str, help='Input file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    structure_type, kwargs = parse_inputs(args.input)
    run_watcon(structure_type, kwargs)


    
