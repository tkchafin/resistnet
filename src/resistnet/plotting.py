"""
Shared plot styling: one palette for every ResistNet figure (matplotlib
and the Plotly report), from a validated colour-blind-safe reference
palette.
"""

# One colour per model family (a validated categorical pair), and neutral
# grey for reference models
COLOURS = {"resistance": "#2a78d6", "samc": "#eb6834",
           "reference": "#8a8986"}
FAMILY_NAMES = {"resistance": "Standard (resistance)", "samc": "SAMC",
                "reference": "Reference model"}

# Ink, grid and axis
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"

# Sequential single-hue ramp (light to dark) for magnitudes. On maps of
# lines, start no lighter than step 250 so low values stay visible.
BLUES = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
         "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
         "#0d366b"]
BLUES_FOR_LINES = BLUES[3:]

FONT = "IBM Plex Sans, Helvetica Neue, Arial, sans-serif"


def plotly_colorscale(colours):
    """A Plotly colour scale through the given colours, evenly spaced."""
    n = len(colours) - 1
    return [[i / n, c] for i, c in enumerate(colours)]


def plotly_layout(**kw):
    """Plotly layout defaults matching the ResistNet palette."""
    axis = dict(gridcolor=GRID, linecolor=AXIS, zeroline=False,
                tickfont=dict(color=TEXT2, size=11),
                title=dict(font=dict(color=TEXT2, size=12)))
    layout = dict(
        font=dict(family=FONT, color=TEXT, size=12),
        paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=60, r=24, t=48, b=48),
        hoverlabel=dict(bgcolor=TEXT, font=dict(color="white", family=FONT,
                                                size=12)),
        legend=dict(font=dict(color=TEXT2, size=11), orientation="h",
                    yanchor="bottom", y=1.02, xanchor="right", x=1),
        xaxis=axis, yaxis=dict(axis),
        title=dict(font=dict(size=14, color=TEXT), x=0, xanchor="left"),
    )
    layout.update(kw)
    return layout
