#!/bin/python

import pandas as pd
import os
import plotly.graph_objects as go

grid_color = '#ddd9d9'
width = 3


# Predefined order of gene names
gene_order = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP','NS']
#FUNCTION TO PLOT COVERAGE PLOTS
def plot_line_coverage_depth(file_data,irma_dir, reference_folder): #file_data is a dictionary where key is foldername, for the data.
    for foldername, cov_data in file_data.items(): #every folder data frame
        print(f'processing {foldername}')
        positions = []
        max_x_position = 0
        fig = go.Figure()
        cov_data.sort(key=lambda x: gene_order.index(x[1]))
        for df, name, color in cov_data:
            # print(df, name)
            df['Position'] += max_x_position

            fig.add_trace(go.Scatter(
                x=df.Position, y=df['Coverage Depth'], name=name,fill = 'tozeroy',
                line=dict(color=color, width=width), mode='lines+text') )
            positions.append(max_x_position)
            # # Update the maximum x position for the next file
            max_x_position = df['Position'].max() + 1  # Add a small gap for separation

            fig.update_layout(
                title = dict(text = f'{foldername} Coverage',x=0.5,xanchor='center', font=dict(family='Arial Black')),

                plot_bgcolor='white',height=600,width=1400,font=dict(size=20, color='black'),
                legend=dict(x=1, y=1, bordercolor='black', borderwidth=1, font=dict(size=18)),
                xaxis=dict(
                    tickvals=[(positions[i] + positions[i+1]) / 2 for i in range(len(positions) - 1)] + [positions[-1]],
                    ticktext=[name for _, name, _ in cov_data]
            ))

            # Update axes
            fig.update_xaxes(
                linecolor='black',tickfont=dict(size=25),tickcolor='black',title='Gene position')
            
                    # Add vertical lines for visual emphasis
            for pos in positions:
                fig.add_vline(x=pos, line=dict(color='gray', width=0.5))

            fig.update_yaxes(
                linecolor='black',gridcolor=grid_color,
                tickfont=dict(size=25),tickcolor='black',title='Read Depth'
                )
            

        fig.write_image(f'{irma_dir}/{reference_folder}/{foldername}_depth.pdf')
        print(f'Completed processing {foldername}')
    #return fig
        