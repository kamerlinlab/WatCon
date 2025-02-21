import os, sys
import argparse

def parse_inputs(filename):
    """
    Parse inputs from input file

    Returns:
    Structure type ('static'/'dynamic'), kwargs dictionary
    """

    #Initialize blank kwargs dict
    kwargs = {}

    analysis_conditions = {}

    with open(filename, 'r') as FILE:
        lines = FILE.readlines()

    for i, line in enumerate(lines):
        if 'Property calculation' in line:
            analysis_conditions = {}
            for analysis_condition in lines[i+1:]:
                if len(analysis_condition) < 1 or 'analysis_selection' in analysis_condition:
                    break
                else:
                    kw = analysis_condition.split(':')[0]
                    kw_value = analysis_condition.split(':')[1].split()[0]
                    analysis_conditions[kw] = kw_value

        #Skip past comments    
        if line.startswith(';') or ':' not in line:
            continue

        #Choose dynamic or static
        elif 'structure_type' in line:
            structure_type = line.split(':')[1].split()[0]

        elif 'make_fastas' in line:
            make_fastas = line.split(':')[1].split()[0]

        #Assign all other kwargs
        else:
            kw = line.split(':')[0]
            kw_value = ' '.join(line.split(':')[1].split()[0].split("#"))
            if kw_value == 'on':
                kw_value = True
            elif kw_value == 'off':
                kw_value = False

            if kw not in analysis_conditions.keys():
                kwargs[kw] = kw_value

    if len(analysis_conditions.keys()) == 0:
        analysis_conditions = 'all'
    kwargs['analysis_conditions'] = analysis_conditions

    if 'water_reference_resids' in kwargs.keys():
        final_resids = []
        resids = kwargs['water_reference_resids'].split(',')
        for resid in resids:
            final_resids.append(int(resid))
        kwargs['water_reference_resids'] = final_resids

    #Change strings to floats/ints
    if 'max_distance' in kwargs.keys():
        kwargs['max_distance'] = float(kwargs['max_distance'])
    
    if 'angle_criteria' in kwargs.keys():
        if kwargs['angle_criteria'] == 'None':
            kwargs['angle_criteria'] = None
        else:
            kwargs['angle_criteria'] = float(kwargs['angle_criteria'])

    if 'active_site_radius' in kwargs.keys():
        kwargs['active_site_radius'] = float(kwargs['active_site_radius'])

    if 'multi_model_pdb' in kwargs.keys():
        if kwargs['multi_model_pdb'] == 'False':
            kwargs['multi_model_pdb'] = False
        else:
            kwargs['multi_model_pdb'] = True
     
    if 'min_cluster_samples' in kwargs.keys():
        kwargs['min_cluster_samples'] = int(kwargs['min_cluster_samples'])

    if 'eps' in kwargs.keys():
        kwargs['eps'] = float(kwargs['eps'])

    if 'num_workers' in kwargs.keys():
        kwargs['num_workers'] = int(kwargs['num_workers'])
    
    if 'water_name' in kwargs.keys():
        if kwargs['water_name'] == 'default':
            kwargs['water_name'] = None
    
    if 'trajectory_name' in kwargs.keys() and kwargs['trajectory_name'] == 'None':
        kwargs['trajectory_name'] = None

    if 'topology_name' in kwargs.keys() and kwargs['topology_name'] == 'None':
        kwargs['topology_name'] = None
    
    return structure_type, kwargs

#NOTE NEED TO MAKE IT SO THAT YOU CAN TELL IT TO MAKE THE FASTAS FROM THE PDBS

def check_conditions(kwargs):
    print('Checking conditions from input file')
    #Consider making some text to check conditions
    pass
    
def run_watcon(structure_type, kwargs):
    if structure_type == 'static':
        from WatCon.generate_static_networks import initialize_network
    else:
        from WatCon.generate_dynamic_networks import initialize_network

    results = initialize_network(**kwargs)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform analysis using WatCon')
    parser.add_argument('--input', type=str, help='Input file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    structure_type, kwargs = parse_inputs(args.input)
    run_watcon(structure_type, kwargs)


    
