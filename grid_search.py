import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_script(pweight, nstart, ncand, gamma):
    command = [
        './scripts/runResistnet.py',
        '-g', '/Users/tyler/projects/resistnet_validation/simulation/output/tc1_s2_30_1.ResistanceMatrix.tsv',
        '-n', '/Users/tyler/projects/resistnet_validation/simulation/template_networks/tc1.net',
        '-c', '/Users/tyler/projects/resistnet_validation/simulation/output/tc1_s2_30_1.coords',
        '-V', '/Users/tyler/projects/resistnet_validation/simulation/selected_vars.txt',
        '-F', '50',
        '-i', '100',
        '-t', '1',
        '--reps', str(10),
        '-P', str(pweight),
        '-G', str(gamma),
        '-S', str(nstart),
        '-C', str(ncand)
    ]
    
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {command}")
        print(result.stderr)
        return None
    return result.stdout

def extract_losses(output):
    pattern = r"Worker \d+: Best Loss = ([\d\.]+)"
    matches = re.findall(pattern, output)
    losses = [float(match) for match in matches]
    return losses

def grid_search(pweight_values, nstart_values, ncand_values, gamma_values):
    results = []

    with ThreadPoolExecutor(max_workers=6) as executor:
        future_to_params = {
            executor.submit(run_script, pweight, nstart, ncand, gamma): (pweight, nstart, ncand, gamma)
            for pweight in pweight_values
            for nstart in nstart_values
            for ncand in ncand_values
            for gamma in gamma_values
        }

        for future in as_completed(future_to_params):
            params = future_to_params[future]
            try:
                output = future.result()
                if output:
                    losses = extract_losses(output)
                    if losses:
                        loss_spread = np.std(losses)
                        results.append((params[0], params[1], params[2], params[3], loss_spread, np.mean(losses)))
                    else:
                        print(f"No losses found for parameters: {params}")
            except Exception as exc:
                print(f"Generated an exception: {params} - {exc}")

    return results

def main():
    # Parameter grid for TPE
    pweight_values = [0.3, 0.5, 0.7, 0.9]
    nstart_values = [10, 30, 50]
    ncand_values = [20, 80, 120]
    gamma_values = [0.1, 0.2, 0.3, 0.4, 0.5]

    results = grid_search(pweight_values, nstart_values, ncand_values, gamma_values)

    if not results:
        print("No valid results found.")
        return

    # Find the best parameter combination
    results = sorted(results, key=lambda x: x[5])  # Sort by mean loss

    # Prepare data for plotting
    pw, ns, nc, g, spreads, means = zip(*results)

    # Data for seaborn
    import pandas as pd
    df = pd.DataFrame({
        'Prior Weight (pweight)': pw,
        'Initial Random Evaluations (nstart)': ns,
        'EI Candidate Points (ncand)': nc,
        'Exploration Factor (gamma)': g,
        'Mean Best Loss': means,
        'Loss Spread': spreads
    })

    # Pairplot to show relationships between all parameters
    sns.pairplot(df, hue='Mean Best Loss', palette='viridis')
    plt.show()

    # Parallel coordinates plot
    from pandas.plotting import parallel_coordinates
    plt.figure(figsize=(12, 6))
    parallel_coordinates(df[['Prior Weight (pweight)', 'Initial Random Evaluations (nstart)', 'EI Candidate Points (ncand)', 'Exploration Factor (gamma)', 'Mean Best Loss']], class_column='Mean Best Loss', colormap=plt.get_cmap("viridis"))
    plt.title('Parallel Coordinates Plot')
    plt.show()

    print("Best parameter combination (pweight, nstart, ncand, gamma):", results[0][:4])
    print("With mean best loss:", results[0][5])

if __name__ == "__main__":
    main()