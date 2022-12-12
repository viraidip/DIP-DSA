motif_search_tab <- tabItem(tabName="motif_search",
  h1("Motif search"),
  fluidRow(
    box(
      title="Input motif",
      width=12,
      textInput(
        inputId="motif",
        label="Type a motif to search for"
      ),
    ),
    box(
      title="Motif matches on sequence",
      width=12,
      "Showing the matches to the motif on the full length sequence",
      plotlyOutput("motif_on_sequence")
    ),
    box(
      title="Table of matches",
      width=12,
      "Showing the matches to the motif as a table",
      dataTableOutput("motif_table")
    ),
  )
)
