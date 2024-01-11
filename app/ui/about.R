about_tab <- tabItem(tabName="about",
  h1("About this application"),
  fluidRow(
    box(
      title="What this application is meant for",
      width=12,
      "This application allows to investigate datasets of deletion-containing",
      "viral genomes (DelVGs). They are defined by a big deletion in the",
      "middle of the original/wild-type sequence. This deletion is described",
      "in the corresponding datasets by the start and the end point of the",
      "deletion."
    ),
    box(
      title="Overview about the different tabs",
      width=12,
      "The application consists of multiple tabs, that provide different",
      "analyses.",
      tags$ol(
        tags$li(tags$b("Add new dataset:"),
          "Allows to upload a new custom dataset. If a strain was used that",
          "is not included right now it has to be added first. For that the",
          "FASTA files of the single segments need to be provided."
        ),
        tags$li(tags$b("Single dataset:"),
          "Allows the user to investigate a single dataset. The read support",
          "cutoff (RSC) can be indivudally set. Additionally, the user ",
          "decides if the NGS counts should be included. "
        ),
        tags$li(tags$b("Multiple datasets:"),
          "Provides similar analyses as the 'single dataset' tab. But in this",
          "case multiple datasets can be compared. "
        ),
        tags$li(tags$b("Dataset intersection:"),
          "In this tab the user can search for DVGs that are present in",
          "multiple datasets. It shows the candidates and the overlap between",
          "the datasets. Additionally, a plot of the NGS counts is given",
          "where the DVGs with the highest occurrence are marked."
        )
      )
    ),
    box(
      title="How to enter a custom dataset?",
      width=12,
      strong("General info:"),
      tags$ol(
        tags$li("Custom datasets need to be in *.csv format."),
        tags$li("They have exactly four columns."),
        tags$li("Ordering of the four columns is crucial!"),
      ),
      strong("Example dataset:"),
      br(),
      "Including correct order and naming of the headers",
      tags$table(border=2, width="100%",
        tags$tbody(
          tags$tr(
            tags$td(""),
            tags$td(strong("Column 1")),
            tags$td(strong("Column 2")),
            tags$td(strong("Column 3")),
            tags$td(strong("Column 4"))
          ),
          tags$tr(
            tags$td(strong("header names")),
            tags$td(tags$i("Segment")),
            tags$td(tags$i("Start")),
            tags$td(tags$i("End")),
            tags$td(tags$i("NGS_read_count"))
          ),
          tags$tr(
            tags$td(strong("column data type")),
            tags$td("character (string)"),
            tags$td("integer"),
            tags$td("integer"),
            tags$td("integer")
          ),
          tags$tr(
            tags$td(strong("description")),
            tags$td("Name of the segment that the DelVG is coming from"),
            tags$td("Start position of the deletion site"),
            tags$td("End position of the deletion site"),
            tags$td("Number of counts in the NGS data of this specific DelVG",
              "candidate"
            )
          )
        )
      )
    ),
    box(
      title="What are 'direct repeats'?",
      width=12,
      "The sequence of the RNA before the starting point and end point can be",
      "the same (is repeated). If this is the case it is called a 'direct",
      "repat'. Direct repeats can be of different length and are disscussed",
      "to be a driving factor in the generation of DelVGs.",
      br(),
      tags$img(src="direct_repeats.png"),
      br(),
      "The actual start and end point of a DIP that has a direct repeat",
      "longer than 0 can not be determined certainly. In the given image a",
      "direct repeat of length 3 is given. There are four possible ways on",
      "how it could have been created during the experiment. If n is the",
      "length of the direct repeat, there are always n+1 options on how it",
      "could have been created."
    ),
    box(
      title="Contact information",
      width=12,
      "To get into contact with the developement team open a new",
      tags$a(href="https://github.com/LohmannJens/DIP-DSA/issues", "issue"),
      "on GitHub. There you can get help with more detailed questions and",
      "come up with new ideas and features."
    )
  )
)
