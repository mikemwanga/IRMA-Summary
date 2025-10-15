#!/bin/python3

from config import *

#PLOT FUNCTIONS

def plot_coverage_depth(data, header,bartitle,max_value,color_scale):
    '''Function to plot coverage by depth'''
    fig = px.imshow(data,text_auto=True,aspect='auto',color_continuous_scale=color_scale,zmin=0, zmax=max_value)
    fig.update_layout(margin = margin,title = header,font_color=font_color,
                      plot_bgcolor = "white",#height=300,
                      coloraxis_colorbar=dict(title=bartitle,orientation='v',#y=-0.3,
                                              title_font=dict(size=12),
                                              title_side='top',thickness=13))
    fig.update_yaxes(tickfont=dict(size=14), title = None, linecolor='black',
                     showticklabels=True,ticks='outside',tickcolor='black')
    fig.update_xaxes(tickfont=dict(size=14),title = 'Gene_Name (Ordered from longest to shortest)', linecolor='black',
                     showticklabels=True,ticks='outside',tickcolor='black')
    fig.update_traces(ygap=0.5,xgap=0.5)
    return fig


#*************************************
print('Plotting figures for coverage by depth')

def plot_figures_log(data):
    max_value = data.max(numeric_only=True).max()
    binsize=20
    num_bins = (len(data) + binsize -1) // binsize #// RETURNS A ROUNDED OFF VALUE AND NOT FLOAT
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize,len(data))
        subset_df = data.iloc[start:stop, :]
        name =f'Subset_{i+1}'
        fig = plot_coverage_depth(subset_df, f'Read Coverage - Log Counts_{name}','Log10', max_value,'OrRd')
        # print(fig)
        # fig.write_html(f'{irma_dir}/coverage-depth-log-counts-{name}.html')
        # pio.write_image(fig, f'coverage_length.pdf', engine='kaleido', format=.pdf')
        fig.write_image(f'{irma_dir}/coverage_length_{name}.pdf', engine='kaleido', format='pdf')

def plot_figures_length(data):
    max_value = data.max(numeric_only=True).max()
    binsize=20
    num_bins = (len(data) + binsize -1) // binsize #// RETURNS A ROUNDED OFF VALUE AND NOT FLOAT
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize,len(data))
        subset_df = data.iloc[start:stop, :]
        name =f'Subset-{i+1}'
        fig = plot_coverage_depth(subset_df, f'Gene Coverage - Length_{name}','Coverage (%)',max_value,'GnBu' )
        fig.write_image(f'{irma_dir}/coverage_depth_{name}.pdf', format='pdf')
        # pio.write_image(fig, f'{irma_dir}/coverage_length_{name}.pdf')


#Plot total reads and assembled 
def plot_bar(data,header):
    '''Function to plot reads assembled relative to total number of reads'''
    fig = px.bar(data, y='sample',x = 'assembled', text="total_reads", range_x=(0,100))
    fig.update_layout(margin = margin,title = header,plot_bgcolor = "white", font_color=font_color)
    fig.update_yaxes(title_font=dict(size=14),tickfont=dict(size=14),tickcolor='black',title=None,
                     linecolor='black',showticklabels=True,ticks='outside')
    fig.update_xaxes(linecolor='black',title_font=dict(size=14),tickfont=dict(size=14),tickcolor='black',
                     title = '% Reads Assembled',gridcolor='#ddd9d9',
                     showticklabels=True,ticks='outside')
    fig.update_traces(marker_color = '#367588',width=0.8,
                      textfont_size=12, textangle=0, textposition="inside", cliponaxis=False)
    return fig

def plot_assembled_reads(data):
    data.sort_values('sample', inplace=True)
    binsize=20
    num_bins = (len(data) + binsize -1) // binsize #// RETURNS A ROUNDED OFF VALUE AND NOT FLOAT
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize,len(data))
        subset_df = data.iloc[start:stop, :]
        name =f'Subset-{i+1}'
        fig = plot_bar(subset_df, f'{name}_total-reads_vs_%-assembled' )
        fig.write_image(f'{irma_dir}/assembled_reads_{name}.pdf')


# Predefined order of gene names
gene_order = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'NS','MP']
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
                x=df.Position, y=df['Consensus_Count'], name=name,fill = 'tozeroy',
                line=dict(color=color, width=width), mode='lines+text') )
            positions.append(max_x_position)
            # # Update the maximum x position for the next file
            max_x_position = df['Position'].max() + 1 # Add a small gap for separation
            #print([(positions[i] + positions[i+1]) / 2 for i in range(len(positions) - 1)] + [positions[-1]])
            fig.update_layout(showlegend=False,
                title = dict(text = f'{foldername} Coverage',x=0.5,xanchor='center', font=dict(family='Arial Black')),

                plot_bgcolor='white',height=600,width=1400,font=dict(size=20, color='black'),
                #legend=dict(x=1, y=1, bordercolor='black', borderwidth=1, font=dict(size=18)),
                legend = None,
                xaxis=dict(#tickmode='linear',
                    tickvals=[(positions[i] + positions[i+1])/2 for i in range(len(positions) - 1)] + [positions[-1]],
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