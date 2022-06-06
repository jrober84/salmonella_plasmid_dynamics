import pandas as pd
import plotly.express as px


def read_sample_data(file):
    df = pd.read_csv(file, sep="\t", header=0)
    return df


def create_sunburst(df):
    path = ['drug_class',
            'gene_id', ]

    fig = px.sunburst(
        data_frame=df,
        path=path,
        values='total',
        maxdepth=3,
        color='proportion',
        color_continuous_scale='RdBu',

    )

    fig.update_layout(
        margin=dict(t=1, l=1, r=1, b=1)
    )
    fig.update_traces(textinfo="label+value")

    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'ResistanceGeneSubburst',
            'height': 600,
            'width': 1200,
            'scale': 3  # Multiply title/legend/axis/canvas sizes by this factor
        },
    }
    fig.update_layout(
        font=dict(
            size=14,
        ),
        coloraxis_colorbar_title='Proportion chromosomal',
        coloraxis_colorbar_x=1
        # colorbar={"orientation": "h", "x": 0.5, "yanchor": "middle", "y": 0.1},

    )

    fig.show()
    return fig


def scater_plot(df):
    fig = px.scatter(df[df['plasmid'] > 10], x="plasmid entropy", y="serovar entropy",
                     size="total",
                     color='human proportion',
                     hover_name="gene_id", log_x=False, size_max=60, labels={
            "serovar entropy": "Serovar entropy",
            "plasmid entropy": "Plasmid entropy",
            "human proportion": "Human proportion",

        }, )
    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 600,
            'width': 1200,
            'scale': 3  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        font=dict(
            size=14,
        )
    )
    fig.show(config=config)
    return


def scater_plot_chao1(df):
    fig = px.scatter(df[df['plasmid_chao1'] >= 0], x="total_samples", y="plasmid_chao1",
                     color='proportion_resistant',
                     hover_name="serovar", log_x=False, trendline="ols", )
    fig.update_layout(
        font=dict(
            size=14,
        )
    )
    fig.show()
    return


def scater_plot_mobility(df):
    df = df[df['serovar_entropy'] >= 0]
    fig = px.scatter(df, x="serovar_entropy", y="total_samples",
                     color='overall_mobility',
                     labels={
                         "serovar_entropy": "Serovar entropy",
                         "total_samples": "log(10) Total cluster members",
                         'overall_mobility': "Overall mobility"
                     },
                     hover_name="plasmid_id", log_y=True, size_max=60, marginal_x="violin")
    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 600,
            'width': 1200,
            'scale': 3  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        font=dict(
            size=14,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ))
    fig.show(config=config)
    return


def scater_plot_plasmid_resistance(df):
    fig = px.scatter(df, x="proportion_resistant", y="serovar_entropy",
                     size='total_samples', color='overall_mobility',
                     labels={
                         "proportion_resistant": "Proportion resistant",
                         "serovar_entropy": "Serovar entropy",
                         'overall_mobility': "Overall mobility"
                     },
                     hover_name="plasmid_id", log_y=False, size_max=60)
    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 600,
            'width': 1200,
            'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        font=dict(
            size=14,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1)
    )
    fig.show(config=config)
    return

def serovar_plasmid_frac(df):

    fig = px.bar(df, x="serovar",
                 y=["positive", "negative"], log_y=True,
                 color_discrete_map={
                     "positive": "#2978A0", "negative": "#F17300"
                 },
                 )
    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'serovar_plasmid_frac',
            'height': 600,
            'width': 1200,
            'scale': 2  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        legend_title="Plasmid presence",
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",

    )
    fig.update_yaxes(title="log(10) Sample count")
    fig.update_xaxes(title="Serovar")
    fig.show(config=config)

def collapse_serovars(df,num_serovars=10):
    df = df.sort_values('total_samples',ascending = False)
    counts = {}
    num_labels = 0
    positive = 0
    negative = 0
    for index,row in df.iterrows():
        serovar = row.serovar
        count_plasmid_negative_samples = row.count_plasmid_negative_samples
        count_plasmid_positive_samples = row.count_plasmid_positive_samples
        if index < num_serovars:
            counts[serovar] = {'serovar':serovar,'positive':count_plasmid_positive_samples,'negative':count_plasmid_negative_samples}
        else:
            num_labels+=1
            positive+=count_plasmid_positive_samples
            negative+=count_plasmid_negative_samples
    counts["{}_others".format(num_labels)] = {'serovar':"{}_others".format(num_labels),'positive':positive,'negative':negative}
    return pd.DataFrame.from_dict(counts,orient='index')





plasmid_mob_serovar_df = read_sample_data('/Users/jrobertson/GoogleDrive/Pycharm/PycharmProjects/salmonella_plasmid_dynamics/results/plasmid_serovar_entropy.txt')
scater_plot_mobility(plasmid_mob_serovar_df)
scater_plot_plasmid_resistance(plasmid_mob_serovar_df)
df = read_sample_data(
    '/Users/jrobertson/GoogleDrive/Pycharm/PycharmProjects/salmonella_plasmid_dynamics/results/genes.txt')
create_sunburst(df)
gene_df = read_sample_data(
    '/Users/jrobertson/GoogleDrive/Pycharm/PycharmProjects/salmonella_plasmid_dynamics/results/salmonella_res_genes.txt')
scater_plot(gene_df)

serovar_df = read_sample_data('/Users/jrobertson/GoogleDrive/Pycharm/PycharmProjects/salmonella_plasmid_dynamics/results/serovar.txt')
serovar_abundance = collapse_serovars(serovar_df)

serovar_plasmid_frac(serovar_abundance)
