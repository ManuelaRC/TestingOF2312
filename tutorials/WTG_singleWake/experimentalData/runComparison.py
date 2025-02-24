import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# Step 1: Read from CSV file (simulation data)
file_path = 'openfoam_dataInPlaneZ90m.csv'
data_kefp = pd.read_csv(file_path)

# 
file_path_laan25FP = 'expData_Laan_NREL_U_HighTI_2-5D_FP.csv'
data_laan25FP = pd.read_csv(file_path_laan25FP)

file_path_laan5FP = 'expData_Laan_NREL_U_HighTI_5D_FP.csv'
data_laan5FP = pd.read_csv(file_path_laan5FP)

file_path_laan75FP= 'expData_Laan_NREL_U_HighTI_7-5D_FP.csv'
data_laan75FP = pd.read_csv(file_path_laan75FP)

# 
file_path_laan25LES = 'expData_Laan_NREL_U_HighTI_2-5D_LES.csv'
data_laan25LES = pd.read_csv(file_path_laan25LES)

file_path_laan5LES = 'expData_Laan_NREL_U_HighTI_5D_LES.csv'
data_laan5LES = pd.read_csv(file_path_laan5LES)

file_path_laan75LES= 'expData_Laan_NREL_U_HighTI_7-5D_LES.csv'
data_laan75LES = pd.read_csv(file_path_laan75LES)

# Step 3: Create a list with values
values_list = [2.5, 5, 7.5]

# Step 4: User input for diameter value
diameter = 126

# Step 5: User input for disc location in x and y axes
x_location = 6.30
y_location = 0

font_size = 18

# Find closest match for x and y locations
closest_x_location = data_kefp['Points_0'].iloc[(data_kefp['Points_0'] - x_location).abs().idxmin()]
closest_y_location = data_kefp['Points_1'].iloc[(data_kefp['Points_1'] - y_location).abs().idxmin()]

print(f"Closest x-axis location: {closest_x_location}")
print(f"Closest y-axis location: {closest_y_location}")

fig, axs = plt.subplots(1, 3, figsize=(20, 6))

# Step 6-8: Plot arcs for each item in the values list
for value in values_list:
    # Calculate radius for the arc
    radius = value * diameter

    # Calculate angles for arc (from -30 degrees to +30 degrees)
    angles = np.linspace(-40, 40, 100)

    # Calculate x and y coordinates for the arc
    x_arc = closest_x_location + radius * np.cos(np.radians(angles))
    y_arc = closest_y_location + radius * np.sin(np.radians(angles))

    # Create DataFrames for x_arc and y_arc
    arc_df = pd.DataFrame({'x_arc': x_arc, 'y_arc': y_arc})

    # Find nearest points in the dataset to the arc coordinates
    distances = np.sqrt(((data_kefp[['Points_0', 'Points_1']].values[:, None] - arc_df.values) ** 2).sum(axis=2))
    nearest_points_index = distances.argmin(axis=0)
    nearest_points = data_kefp.iloc[nearest_points_index]

    # Extract 'U_0' values for the nearest points and scale by 1/8
    extracted_U_0 = nearest_points['U_0'].values * (1 / 8)

    # Apply moving average filter to extracted_U_0, leaving boundary points unchanged
    window_size = 10  # Adjust window size as needed
    # Smooth interior points
    smoothed_U_0_interior = np.convolve(extracted_U_0[1:-1], np.ones(window_size)/window_size, mode='same')

    # Concatenate the boundary points with the smoothed interior points
    smoothed_U_0 = np.concatenate(([extracted_U_0[0]], smoothed_U_0_interior, [extracted_U_0[-1]]))

    # Plotting the figure for each value in the list
    
    if value == 2.5:  # Plot experimental data only for value = 2.5
        axs[0].plot(data_laan25LES.iloc[:, 0], data_laan25LES.iloc[:, 1], label='LES 6x10 min', linestyle='-', color='blue')
        axs[0].plot(data_laan25FP.iloc[:, 0], data_laan25FP.iloc[:, 1], label='k-e-fp M.P. van der Laan', linestyle='-', color='green')
        axs[0].plot(angles, smoothed_U_0, label=f"k-e-fp OpenFOAM", linestyle='-', color='black')
        # Set x-axis and y-axis limits
        axs[0].set_xlim(-30, 30)
        axs[0].set_ylim(0.5, 1.1)
        axs[0].set_title("2.5 D", fontstyle = "italic", fontsize=font_size)
        axs[0].set_ylabel("$U/U_{\mathrm{ref}}$ [-]", fontsize=font_size)
        axs[0].tick_params(axis='both', which='major', labelsize=font_size-2)

    if value == 5:  # Plot experimental data only for value = 5
        axs[1].plot(data_laan5LES.iloc[:, 0], data_laan5LES.iloc[:, 1], label='LES 6x10 min', linestyle='-', color='blue')
        axs[1].plot(data_laan5FP.iloc[:, 0], data_laan5FP.iloc[:, 1], label='k-e-fp M.P. van der Laan', linestyle='-', color='green')
        axs[1].plot(angles, smoothed_U_0, label=f"k-e-fp OpenFOAM", linestyle='-', color='black')
        # Set x-axis and y-axis limits
        axs[1].set_xlim(-30, 30)
        axs[1].set_ylim(0.5, 1.1)
        axs[1].set_title("5 D", fontstyle = "italic", fontsize=font_size)
        axs[1].tick_params(axis='both', which='major', labelsize=font_size-2)

    if value == 7.5:  # Plot experimental data only for value = 7.5
        axs[2].plot(data_laan75FP.iloc[:, 0], data_laan75FP.iloc[:, 1], label='k-e-fp M.P. van der Laan', linestyle='-', color='green')
        axs[2].plot(data_laan75LES.iloc[:, 0], data_laan75LES.iloc[:, 1], label='LES 6x10 min', linestyle='-', color='blue')
        axs[2].plot(angles, smoothed_U_0, label=f"k-e-fp OpenFOAM", linestyle='-', color='black')

        # Set x-axis and y-axis limits
        axs[2].set_xlim(-20, 20)
        axs[2].set_ylim(0.5, 1.1)
        axs[2].set_title("7.5 D", fontstyle = "italic", fontsize=font_size)
        axs[2].tick_params(axis='both', which='major', labelsize=font_size-2)
    # Add legend
    axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), shadow=False, ncol=4, fontsize=font_size)

# Add grid
for ax in axs:
    ax.grid(True, linestyle=':', color='gray')    

# Show all the figures
plt.savefig("ComparisonU.png", dpi=600, bbox_inches='tight')
plt.show()

## ------------------------ TI PLOTS ------------------------- ##

# 
file_path_laan25FP = 'expData_Laan_NREL_TI_HighTI_2-5D_FP.csv'
data_laan25FP = pd.read_csv(file_path_laan25FP)

file_path_laan5FP = 'expData_Laan_NREL_TI_HighTI_5D_FP.csv'
data_laan5FP = pd.read_csv(file_path_laan5FP)

file_path_laan75FP= 'expData_Laan_NREL_TI_HighTI_7-5D_FP.csv'
data_laan75FP = pd.read_csv(file_path_laan75FP)

# 
file_path_laan25LES = 'expData_Laan_NREL_TI_HighTI_2-5D_LES.csv'
data_laan25LES = pd.read_csv(file_path_laan25LES)

file_path_laan5LES = 'expData_Laan_NREL_TI_HighTI_5D_LES.csv'
data_laan5LES = pd.read_csv(file_path_laan5LES)

file_path_laan75LES= 'expData_Laan_NREL_TI_HighTI_7-5D_LES.csv'
data_laan75LES = pd.read_csv(file_path_laan75LES)

fig, axs = plt.subplots(1, 3, figsize=(20, 5))

# Step 6-8: Plot arcs for each item in the values list
for value in values_list:
    # Calculate radius for the arc
    radius = value * diameter

    # Calculate angles for arc (from -30 degrees to +30 degrees)
    angles = np.linspace(-40, 40, 100)

    # Calculate x and y coordinates for the arc
    x_arc = closest_x_location + radius * np.cos(np.radians(angles))
    y_arc = closest_y_location + radius * np.sin(np.radians(angles))

    # Create DataFrames for x_arc and y_arc
    arc_df = pd.DataFrame({'x_arc': x_arc, 'y_arc': y_arc})

    # Find nearest points in the dataset to the arc coordinates
    distances = np.sqrt(((data_kefp[['Points_0', 'Points_1']].values[:, None] - arc_df.values) ** 2).sum(axis=2))
    nearest_points_index = distances.argmin(axis=0)
    nearest_points = data_kefp.iloc[nearest_points_index]


    # Extract 'U_0' values for the nearest points and scale by 1/8
    extracted_TI = nearest_points['TI'].values

    # Apply moving average filter to extracted_U_0, leaving boundary points unchanged
    window_size = 10  # Adjust window size as needed
    # Smooth interior points
    smoothed_TI_interior = np.convolve(extracted_TI[1:-1], np.ones(window_size)/window_size, mode='same')

    # Concatenate the boundary points with the smoothed interior points
    smoothed_TI = np.concatenate(([extracted_TI[0]], smoothed_TI_interior, [extracted_TI[-1]]))
    # Plotting the figure for each value in the list
    
    if value == 2.5:  # Plot experimental data only for value = 2.5
        axs[0].plot(data_laan25LES.iloc[:, 0], data_laan25LES.iloc[:, 1], label='LES 6x10 min', linestyle='-', color='blue')
        axs[0].plot(data_laan25FP.iloc[:, 0], data_laan25FP.iloc[:, 1], label='k-e-fp M.P. van der Laan', linestyle='-', color='green')
        axs[0].plot(angles, smoothed_TI, label=f"k-e-fp OpenFOAM", linestyle='-', color='black')
        # Set x-axis and y-axis limits
        axs[0].set_xlim(-30, 30)
        axs[0].set_ylim(0, 0.2)
        axs[0].set_ylabel("$TI$ [-]", fontsize=font_size)
        axs[0].tick_params(axis='both', which='major', labelsize=font_size-2)
    if value == 5:  # Plot experimental data only for value = 5
        axs[1].plot(data_laan5LES.iloc[:, 0], data_laan5LES.iloc[:, 1], label='LES 6x10 min', linestyle='-', color='blue')
        axs[1].plot(data_laan5FP.iloc[:, 0], data_laan5FP.iloc[:, 1], label='k-e-fp M.P. van der Laan', linestyle='-', color='green')
        axs[1].plot(angles, smoothed_TI, label=f"k-e-fp OpenFOAM", linestyle='-', color='black')
        # Set x-axis and y-axis limits
        axs[1].set_xlim(-30, 30)
        axs[1].set_ylim(0, 0.2)
        #axs[1].set_title("5 D", fontstyle = "italic")
        axs[1].set_xlabel("Relative Wind Direction [º]", fontsize=font_size)
        axs[1].tick_params(axis='both', which='major', labelsize=font_size-2)
    if value == 7.5:  # Plot experimental data only for value = 7.5
        axs[2].plot(data_laan75FP.iloc[:, 0], data_laan75FP.iloc[:, 1], label='k-e-fp M.P. van der Laan', linestyle='-', color='green')
        axs[2].plot(data_laan75LES.iloc[:, 0], data_laan75LES.iloc[:, 1], label='LES 6x10 min', linestyle='-', color='blue')
        axs[2].plot(angles, smoothed_TI, label=f"k-e-fp OpenFOAM", linestyle='-', color='black')
        # Set x-axis and y-axis limits
        axs[2].set_xlim(-20, 20)
        axs[2].set_ylim(0, 0.2)
        axs[2].tick_params(axis='both', which='major', labelsize=font_size-2)
        #axs[2].set_title("7.5 D", fontstyle = "italic")


# Add grid
for ax in axs:
    ax.grid(True, linestyle=':', color='gray')    

# Show all the figures
plt.savefig("ComparisonTI.png", dpi=600, bbox_inches='tight')
plt.show()
