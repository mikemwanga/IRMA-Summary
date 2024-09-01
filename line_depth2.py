import pandas as pd
import os
import plotly.graph_objects as go

# Directory and settings
DIR = '/storage/homefs/mm23m894/Assembly_Try/Assembly/IRMA_Default/H17'
grid_color = '#ddd9d9'
width = 3

# File extension to search for
file_extension = 'coverage.txt'

# Predefined order of gene names
gene_order = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP','NS']

# Initialize figure
fig = go.Figure()
 # Read the file

# Choose color based on gene name (if needed)
color_map = {
                'HA': 'blue','NA': '#075d1c','MP': 'black','PA': '#ff21b2','PB2': '#5a189a','PB1': '#9e5231',
                'NP': 'red','NS': '#339dff',
            }
# Initialize the maximum position tracker and lists for positions and names
max_x_position = 0
positions = []  # List to store the starting position of each file's data
file_names = [] # List to store the names of the files for labeling
file_data = []  # List to store dataframes and their colors

# Loop through files in directory
for file in os.listdir(os.path.join(DIR, 'tables')):
    if file.endswith(file_extension):
        # Extract gene name from file name (e.g., 'A_MP-coverage.txt' -> 'MP')
        gene_name = file.split('_')[1].split('-')[0]
        
        # Check if the gene name is in the predefined order
        if gene_name in gene_order:
            name = gene_name  # Use gene name as trace name
            color = None  # Color could be set dynamically or fixed if desired
            
            
            color = color_map.get(name, 'gray')  # Default to gray if not found

            path = os.path.join(DIR, 'tables', file)
            df = pd.read_table(path, usecols=['Reference_Name', 'Position', 'Coverage Depth'])
            
            # Append file data to the list
            file_data.append((df, name, color))

# Sort files according to predefined gene order
file_data.sort(key=lambda x: gene_order.index(x[1]))

# Plot each file after sorting
for df, name, color in file_data:
    # Adjust Position to create horizontal spacing between files
    df['Position'] += max_x_position
    
    # Add trace for this file's data
    fig.add_trace(
        go.Scatter(
            x=df.Position, y=df['Coverage Depth'], name=name,
            line=dict(color=color, width=width), mode='lines+text',
            # text=[name if j == len(df) - 1 else "" for j in range(len(df))],
            # textposition="top right"
        )
    )
    
    # Store the starting position for each file
    positions.append(max_x_position)
    
    # Update the maximum x position for the next file
    max_x_position = df['Position'].max() + 1  # Add a small gap for separation

# Update layout with x-axis labels
fig.update_layout(
    plot_bgcolor='white',height=600,width=1400,
    font=dict(size=20, color='black'),
    legend=dict(x=1, y=1, bordercolor='black', borderwidth=1, font=dict(size=18)),
    xaxis=dict(
        tickvals=[(positions[i] + positions[i+1]) / 2 for i in range(len(positions) - 1)] + [positions[-1]],
        ticktext=[name for _, name, _ in file_data]
    )
)

# Update axes
fig.update_xaxes(
    linecolor='black',tickfont=dict(size=25),tickcolor='black',title='Gene position'
)
fig.update_yaxes(
    linecolor='black',gridcolor=grid_color,
    tickfont=dict(size=25),tickcolor='black',title='Read Depth'
)

# Add vertical lines for visual emphasis
for pos in positions:
    fig.add_vline(x=pos, line=dict(color='gray', width=0.5))

fig.show()
