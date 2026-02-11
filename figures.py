#!/bin/python3

from config import *

# PLOT FUNCTIONS

# MERGE CHANGE: centralize figure aesthetics for consistent output style.
def figure_aesthetics(fig):
    fig.update_xaxes(showticklabels=True, ticks='outside', tickcolor=gridcolor, title=None)
    fig.update_yaxes(showticklabels=True, ticks='outside', tickcolor=gridcolor)
    fig.update_layout(
        margin=margin,font_color=black_color,plot_bgcolor=bgcolor,font_size=textsize,
        legend=dict(orientation='h', y=1.1, x=0.2, title=''))
    return fig


def plot_read_distribution(data, yvalue1, yvalue2, datatype):
    """Plot stacked bars for assembled and non-assembled reads."""
    fig = px.bar(data, y='sample', x=[yvalue1, yvalue2], color_discrete_map=read_colormap)
    fig = figure_aesthetics(fig)
    fig.update_xaxes(title=datatype)
    fig.write_image(f'{figures_dir}/read_distribution-{datatype}.pdf')


def plot_coverage_depth(data, header, bartitle, max_value, color_scale):
    """Function to plot coverage by depth."""
    fig = px.imshow(
        data,
        text_auto=True, aspect='auto',
        color_continuous_scale=color_scale,
        zmin=0,zmax=max_value,
    )
    fig.update_layout(
        margin=margin,title=header, font_color=textcolor, plot_bgcolor=bgcolor,
        coloraxis_colorbar=dict(title=bartitle, orientation='v', thickness=13),
    )
    fig.update_yaxes(
        tickfont=dict(size=textsize),
        title=None,
        linecolor=textcolor,
        showticklabels=True,
        ticks='outside',
        tickcolor=textcolor,
    )
    fig.update_xaxes(
        tickfont=dict(size=textsize),
        title='Gene_Name (Ordered from longest to shortest)',
        linecolor=gridcolor,
        showticklabels=True,
        ticks='outside',
        tickcolor='black',
    )
    fig.update_traces(ygap=0.5, xgap=0.5)
    return fig


print('Plotting figures for coverage by depth')


def plot_figures_log(data):
    max_value = data.max(numeric_only=True).max()
    binsize = 20
    num_bins = (len(data) + binsize - 1) // binsize
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize, len(data))
        subset_df = data.iloc[start:stop, :]
        name = f'Subset_{i+1}'
        fig = plot_coverage_depth(subset_df, f'Read Coverage - Log Counts_{name}', 'Log10', max_value, 'OrRd')
        fig.write_image(f'{figures_dir}/coverage-depth-log-counts-{name}.pdf')


def plot_figures_length(data):
    max_value = data.max(numeric_only=True).max()
    binsize = 20
    num_bins = (len(data) + binsize - 1) // binsize
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize, len(data))
        subset_df = data.iloc[start:stop, :]
        name = f'Subset-{i+1}'
        fig = plot_coverage_depth(subset_df, f'Gene Coverage - Length_{name}', 'Coverage (%)', max_value, 'GnBu')
        fig.write_image(f'{figures_dir}/coverage_length_{name}.pdf')


# Plot total reads and assembled
def plot_bar(data, header):
    """Function to plot reads assembled relative to total number of reads."""
    xcol = 'prop_assembled' if 'prop_assembled' in data.columns else 'assembled'
    fig = px.bar(data, y='sample', x=xcol, text='total_reads', range_x=(0, 100))
    fig.update_layout(margin=margin, title=header, plot_bgcolor=bgcolor, font_color=textcolor)
    fig.update_yaxes(
        tickfont=dict(size=textsize),
        tickcolor=tickfontcolor,
        title=None,
        linecolor=gridcolor,
        showticklabels=True,
        ticks='outside',
    )
    fig.update_xaxes(
        linecolor=gridcolor,
        tickfont=dict(size=textsize),
        tickcolor=tickfontcolor,
        title='% Reads Assembled',
        gridcolor=gridcolor,
        showticklabels=True,
        ticks='outside',
    )
    fig.update_traces(
        marker_color='#367588',
        width=0.8,
        textfont_size=textsize,
        textangle=0,
        textposition='inside',
        cliponaxis=False,
    )
    return fig


def plot_assembled_reads(data):
    data.sort_values('sample', inplace=True)
    binsize = 20
    num_bins = (len(data) + binsize - 1) // binsize
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize, len(data))
        subset_df = data.iloc[start:stop, :]
        name = f'Subset-{i+1}'
        fig = plot_bar(subset_df, f'{name}_total-reads_vs_%-assembled')
        fig.write_image(f'{figures_dir}/assembled_reads_{name}.pdf')


# MERGE CHANGE: safe shared exporter with explicit index control.
def export_data(data, name, index=True):
    data.to_csv(f'{summary_dir}/{name}', sep='\t', index=index)


# Predefined order of gene names
gene_order = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']


def plot_line_coverage_depth(file_data, irma_dir, reference_folder):
    """Plot per-sample concatenated line coverage plot."""
    for foldername, cov_data in file_data.items():
        print(f'processing {foldername}')
        segment_boundaries = []
        segment_centers = []
        max_x_position = 0
        fig = go.Figure()
        cov_data.sort(key=lambda x: gene_order.index(x[1]))

        for df, name, color in cov_data:
            df['Position'] += max_x_position
            segment_start = df['Position'].min()
            fig.add_trace(
                go.Scatter(
                    x=df.Position,
                    y=df['Consensus_Count'],
                    name=name,
                    fill='tozeroy',
                    line=dict(color=color, width=width),
                    mode='lines+text',
                )
            )
            # CHANGE: center tick labels for all segments (including MP/NS) to avoid border-aligned labels.
            segment_end = df['Position'].max()
            segment_centers.append((segment_start + segment_end) / 2)
            max_x_position = df['Position'].max() + 1
            segment_boundaries.append(max_x_position - 0.5)
            fig.update_layout(
                showlegend=False,
                title=dict(text=f'{foldername} Coverage', x=0.5, xanchor='center', font=dict(family=fontfamily)),
                plot_bgcolor=bgcolor,
                height=600,
                width=1400,
                font=dict(size=textsize, color=textcolor),
                legend=None,
                xaxis=dict(
                    tickvals=segment_centers,
                    ticktext=[segment for _, segment, _ in cov_data],
                ),
            )
            fig.update_xaxes(linecolor='black', tickfont=dict(size=textsize), tickcolor='black', title='Gene position')
            for pos in segment_boundaries[:-1]:
                fig.add_vline(x=pos, line=dict(color='gray', width=0.5))
            fig.update_yaxes(
                linecolor='black',
                gridcolor=gridcolor,
                tickfont=dict(size=textsize),
                tickcolor='black',
                title='Read Depth',
            )

        # CHANGE: save individual sample depth plots in dedicated subfolder.
        fig.write_image(f'{individual_depth_dir}/{foldername}_depth.pdf')
        print(f'Completed processing {foldername}')


# MERGE CHANGE: new faceted all-samples plot.
def plot_all_line_plots(data):
    """Plot coverage line graphs for all samples and genes."""
    if data.empty:
        print('No data for coverage-lines-log-counts plot. Skipping.')
        return

    order = ['PB2', 'PB1', 'PA', 'HA', 'NP', '-NA', 'MP', 'NS']
    data = data.copy()
    data.rename(columns={'Reference_Name': 'gene', 'sample_name': 'sample'}, inplace=True)

    fig = px.line(
        data,
        x='Position',
        y='depth_log',
        facet_col='gene',
        color='sample',
        category_orders={'gene': order},
        color_discrete_sequence=px.colors.qualitative.Dark24,
        facet_col_wrap=4,
        range_y=[0, 5],
        # CHANGE: increase facet row spacing to reduce crowding.
        facet_row_spacing=0.28,
        facet_col_spacing=0.04,
    )

    fig = figure_aesthetics(fig)
    fig.update_yaxes(linecolor=gridcolor, mirror=True)
    fig.update_xaxes(linecolor=gridcolor, mirror=True, matches=None)
    # CHANGE: move legend to right side, vertical layout, and smaller legend font.
    fig.update_layout(
        font_size=textsize,
        margin=dict(l=40, r=220, t=50, b=20),
        legend=dict(
            orientation='v',
            y=1,
            yanchor='top',
            x=1.02,
            xanchor='left',
            font=dict(size=max(8, textsize - 6)),
        ),
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split('=')[-1]))
    fig.write_image(f'{figures_dir}/coverage-lines-log-counts.pdf')
