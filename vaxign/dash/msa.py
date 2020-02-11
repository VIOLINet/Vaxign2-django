import dash
import dash_bio as dashbio
import dash_core_components as dcc
import dash_html_components as html

from django_plotly_dash import DjangoDash

app = DjangoDash('msa')

app.layout = html.Div(
    id = 'alignment-viewer-container',
    children = [
        html.Div(id = 'sequence_id', style = {"display":"none"}),
        html.Div(id = 'alignment-viewer'),
    ]
)

@app.callback(dash.dependencies.Output('alignment-viewer', 'children'),
              [dash.dependencies.Input('sequence_id', 'children')])
def display_msa(msa):
    
    return dashbio.AlignmentChart(
        id='alignment-viewer',
        data=open(msa).read(),
        extension='clustal',
    )

if __name__ == '__main__':
    app.run_server(debug=True)