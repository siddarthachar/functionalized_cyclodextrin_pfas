import subprocess
import matplotlib.pyplot as plt

# Function to extract data from .edr file using gmx energy
def extract_energy_data(edr_file, output_file, properties):
    # Convert properties to GROMACS format (separate by newlines)
    properties_input = "\n".join(map(str, properties)) + "\n"
    
    # Run gmx energy
    with open(output_file, 'w') as out:
        subprocess.run(
            ["gmx", "energy", "-f", edr_file, "-o", output_file],
            input=properties_input.encode(),
            stdout=out,
            stderr=subprocess.PIPE
        )

# Specify your .edr file and output .xvg file names
edr_file = "sys-emin.edr"
output_file = "energy_emin_data.xvg"

# Select properties to extract (1=Potential, 2=Kinetic, 3=Temperature, etc.)
properties = [1, 2, 3, 6]  # Adjust based on properties available in .edr file

# Extract data
extract_energy_data(edr_file, output_file, properties)

# Parse the output file to retrieve energy, temperature, pressure data
time = []
potential_energy = []
kinetic_energy = []
temperature = []
pressure = []

with open(output_file, 'r') as f:
    for line in f:
        if not line.startswith(('@', '#')):
            data = line.split()
            time.append(float(data[0]))
            potential_energy.append(float(data[1]))
            kinetic_energy.append(float(data[2]))
            temperature.append(float(data[3]))
            pressure.append(float(data[4]))

# Plot the results
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot Potential Energy
axes[0, 0].plot(time, potential_energy, label="Potential Energy")
axes[0, 0].set_xlabel("Time (ps)")
axes[0, 0].set_ylabel("Potential Energy (kJ/mol)")
axes[0, 0].legend()

# Plot Kinetic Energy
axes[0, 1].plot(time, kinetic_energy, label="Kinetic Energy")
axes[0, 1].set_xlabel("Time (ps)")
axes[0, 1].set_ylabel("Kinetic Energy (kJ/mol)")
axes[0, 1].legend()

# Plot Temperature
axes[1, 0].plot(time, temperature, label="Temperature")
axes[1, 0].set_xlabel("Time (ps)")
axes[1, 0].set_ylabel("Temperature (K)")
axes[1, 0].legend()

# Plot Pressure
axes[1, 1].plot(time, pressure, label="Pressure")
axes[1, 1].set_xlabel("Time (ps)")
axes[1, 1].set_ylabel("Pressure (bar)")
axes[1, 1].legend()

plt.tight_layout()
plt.show()
