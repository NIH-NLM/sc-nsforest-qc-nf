import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import typer
from . import prep_medians, prep_binary_scores, run_nsforest

app = typer.Typer()
app.command(name="prep-medians")(prep_medians.run)
app.command(name="prep-binary-scores")(prep_binary_scores.run)
app.command(name="run-nsforest")(run_nsforest.run)

if __name__ == "__main__":
    app()

