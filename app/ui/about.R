about_tab <- tabItem(tabName="about",
  h1("About this application"),
  fluidRow(
    box(
      title="What this application is meant for",
      width=12,
      "This application allows to investigate datasets of deletion-containing",
      "viral genomes (DelVGs) for influenza viruses. Influenza DelVGs contain",
      "large, usually single internal deletions in their vRNA sequence. This",
      "deletion is defined in the corresponding datasets by the start and the",
      "end point of the deletion. Start corresponds to the 5'-end and End to",
      "the 3'-end of the vRNA sequence."
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
          "cutoff (RSC), which is the minimum NGS count for a DelVG to be ",
          "included, can be set individually. In addition, the user decides",
          "if the DelVGs are weighted by their NGS count (unflattened) or not",
          "(flattened)",
        ),
        tags$li(tags$b("Multiple datasets:"),
          "Provides the same analyses as the 'single dataset' tab. But in",
          "this case multiple datasets can be compared. "
        ),
        tags$li(tags$b("Dataset intersection:"),
          "In this tab the user can search for DelVGs that are present in",
          "multiple datasets. It shows the candidates and the overlap between",
          "the datasets. Additionally, a plot of the NGS counts is given",
          "where the DelVGs with the highest occurrence are marked."
        )
      )
    ),
    box(
      title="Verfied datasets",
      width=12,
      "The following datasets were added by us to the applicaten. All others",
      "are uploaded by other users and we cannot gurantee for their quality.",

      tags$table(border=2, width="100%",
        tags$tbody(
          tags$tr(
            tags$td(strong("Dateset name")),
            tags$td(strong("Strain")),
            tags$td(strong("Publication"))
          ),
          tags$tr(
            tags$td("Alnaji2021"),
            tags$td("A/Puerto Rico/8/1934"),
            tags$td(tags$a(href="https://doi.org/10.1128/mBio.02959-21",
              target="_blank",
              "'Influenza A Virus Defective Viral Genomes Are Inefficiently
              Packaged into Virions Relative to Wild-Type Genomic RNAs.'")
            )
          ),
          tags$tr(
            tags$td("Pelz2021"),
            tags$td("A/Puerto Rico/8/1934"),
            tags$td(tags$a(href="https://doi.org/10.1128/JVI.01174-21",
              target="_blank",
              "'Semi-continuous Propagation of Influenza A Virus and Its
              Defective Interfering Particles: Analyzing the Dynamic
              Competition To Select Candidates for Antiviral Therapy.'")
            )
          ),
          tags$tr(
            tags$td("Wang2023"),
            tags$td("A/Puerto Rico/8/1934"),
            tags$td(tags$a(href="https://doi.org/10.1128/jvi.00493-23",
              target="_blank",
              "'Influenza Defective Interfering Virus Promotes Multiciliated
              Cell Differentiation and Reduces the Inflammatory Response in
              Mice.'")
            )
          ),
          tags$tr(
            tags$td("Wang2020"),
            tags$td("A/Puerto Rico/8/1934"),
            tags$td(tags$a(href="https://doi.org/10.1128/mbio.02880-19",
              target="_blank",
              "'Cell-to-Cell Variation in Defective Virus Expression and
              Effects on Host Responses during Influenza Virus Infection'")
            )
          ),
          tags$tr(
            tags$td("Zhuravlev2020"),
            tags$td("A/Puerto Rico/8/1934"),
            tags$td(tags$a(href="https://doi.org/10.1016/j.dib.2020.106604",
              target="_blank",
              "'RNA-Seq transcriptome data of human cells infected with
              influenza A/Puerto Rico/8/1934 (H1N1) virus'")
            )
          ),
          tags$tr(
            tags$td("Kupke2020"),
            tags$td("A/Puerto Rico/8/1934"),
            tags$td(tags$a(href="https://doi.org/10.3390/v12010071",
              target="_blank",
              "'Single-Cell Analysis Uncovers a Vast Diversity in Intracellular
              Viral Defective Interfering RNA Content Affecting the Large
              Cell-to-Cell Heterogeneity in  Influenza A Virus Replication.'")
            )
          ),
          tags$tr(
            tags$td("VdHoecke2015"),
            tags$td("A/Puerto Rico/8/1934"),
            tags$td(tags$a(href="https://doi.org/10.1186/s12864-015-1284-z",
              target="_blank",
              "'Analysis of the genetic diversity of influenza A viruses using
              next-generation DNA sequencing.'")
            )
          ),
          tags$tr(
            tags$td("Alnaji2019_Cal07"),
            tags$td("A/California/07/2009"),
            tags$td(tags$a(href="https://doi.org/10.1128/JVI.00354-19",
              target="_blank",
              "'Sequencing Framework for the Sensitive Detection and Precise
              Mapping of Defective Interfering Particle-Associated Deletions
              across Influenza A and B Viruses.'")
            )
          ),
          tags$tr(
            tags$td("Alnaji2019_NC"),
            tags$td("A/New Caledonia/20-JY2/1999"),
            tags$td(tags$a(href="https://doi.org/10.1128/JVI.00354-19",
              target="_blank",
              "'Sequencing Framework for the Sensitive Detection and Precise
              Mapping of Defective Interfering Particle-Associated Deletions
              across Influenza A and B Viruses.'")
            )
          ),
          tags$tr(
            tags$td("Mendes2021"),
            tags$td("A/WSN/1933"),
            tags$td(tags$a(href="https://doi.org/10.1371/journal.ppat.1010125",
              target="_blank",
              "'Library-based analysis reveals segment and length dependent
              characteristics of defective influenza genomes.'")
            )
          ),
          tags$tr(
            tags$td("Boussier2020"),
            tags$td("A/WSN/1933"),
            tags$td(tags$a(href="https://doi.org/10.1261/rna.077529.120",
              target="_blank",
              "'RNA-seq accuracy and reproducibility for the mapping and
              quantification of influenza defective viral genomes.'")
            )
          ),
          tags$tr(
            tags$td("Alnaji2019_Perth"),
            tags$td("A/Perth/16/2009"),
            tags$td(tags$a(href="https://doi.org/10.1128/JVI.00354-19",
              target="_blank",
              "'Sequencing Framework for the Sensitive Detection and Precise
              Mapping of Defective Interfering Particle-Associated Deletions
              across Influenza A and B Viruses.'")
            )
          ),
          tags$tr(
            tags$td("Berry2021_A"),
            tags$td("A/Connecticut/Flu122/2013"),
            tags$td(tags$a(href="https://doi.org/10.1101/2021.07.01.450528",
              target="_blank",
              "'High confidence identification of intra-host single nucleotide
              variants for person-to-person influenza transmission tracking in
              congregate settings'")
            )
          ),
          tags$tr(
            tags$td("Penn2022"),
            tags$td("A/turkey/Turkey/1/2005"),
            tags$td(tags$a(href="https://doi.org/10.1128/jvi.01178-22",
              target="_blank",
              "'Levels of Influenza A Virus Defective Viral Genomes Determine
              Pathogenesis in the BALB/c Mouse Model.'")
            )
          ),
          tags$tr(
            tags$td("Lui2019"),
            tags$td("A/Anhui/1/2013"),
            tags$td(tags$a(href="http://doi.org/10.1080/22221751.2019.1611346",
              target="_blank",
              "'SMRT sequencing revealed the diversity and characteristics of
              defective interfering RNAs in influenza A (H7N9) virus infection.
              '")
            )
          ),
          tags$tr(
            tags$td("Alnaji2019_BLEE"),
            tags$td("B/Lee/1940"),
            tags$td(tags$a(href="https://doi.org/10.1128/JVI.00354-19",
              target="_blank",
              "'Sequencing Framework for the Sensitive Detection and Precise
              Mapping of Defective Interfering Particle-Associated Deletions
              across Influenza A and B Viruses.'")
            )
          ),
          tags$tr(
            tags$td("Berry2021_B"),
            tags$td("B/Victoria/504/2000"),
            tags$td(tags$a(href="https://doi.org/10.1101/2021.07.01.450528",
              target="_blank",
              "'High confidence identification of intra-host single nucleotide
              variants for person-to-person influenza transmission tracking in
              congregate settings'")
            )
          ),
          tags$tr(
            tags$td("Valesano2020_Vic"),
            tags$td("B/Victoria/504/2000"),
            tags$td(tags$a(href="https://doi.org/10.1128/JVI.01710-19",
              target="_blank",
              "'Influenza B Viruses Exhibit Lower Within-Host Diversity than
              Influenza A Viruses in Human Hosts.'")
            )
          ),
          tags$tr(
            tags$td("Sheng2018"),
            tags$td("B/Brisbane/60/2008"),
            tags$td(tags$a(href="https://doi.org/10.1099/jgv.0.001018",
              target="_blank",
              "'Identification and characterization of viral defective RNA
              genomes in influenza B virus.'")
            )
          ),
          tags$tr(
            tags$td("Berry2021_B_Yam"),
            tags$td("B/Yamagata/16/1988"),
            tags$td(tags$a(href="https://doi.org/10.1101/2021.07.01.450528",
              target="_blank",
              "'High confidence identification of intra-host single nucleotide
              variants for person-to-person influenza transmission tracking in
              congregate settings'")
            )
          ),
          tags$tr(
            tags$td("Southgate2019"),
            tags$td("B/Yamagata/16/1988"),
            tags$td(tags$a(href="http://doi.org/10.1093/bioinformatics/btz814",
              target="_blank",
              "'Influenza classification from short reads with VAPOR
              facilitates robust mapping pipelines and zoonotic strain
              detection for routine surveillance applications'")
            )
          ),
          tags$tr(
            tags$td("Valesano2020_Yam"),
            tags$td("B/Yamagata/16/1988"),
            tags$td(tags$a(href="https://doi.org/10.1128/JVI.01710-19",
              target="_blank",
              "'Influenza B Viruses Exhibit Lower Within-Host Diversity than
              Influenza A Viruses in Human Hosts.'")
            )
          )
        )
      ),
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
      "Including correct order and naming of the headers. The header names",
      "start with a capitalised letter.",
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
            tags$td("Name of the segment where the DelVG is originating from"),
            tags$td("Start position of the deletion site"),
            tags$td("End position of the deletion site"),
            tags$td("Number of counts in the NGS data of this specific DelVG",
              "candidate"
            )
          )
        )
      ),
      br(),
      strong("Examplary input file with five DelVGs:"),
      tags$table(border=2, width="100%",
        tags$tbody(
          tags$tr(
            tags$td(strong("Segment")),
            tags$td(strong("Start")),
            tags$td(strong("End")),
            tags$td(strong("NGS_read_count"))
          ),
          tags$tr(
            tags$td("PB2"),
            tags$td(163),
            tags$td(2139),
            tags$td(42)
          ),
          tags$tr(
            tags$td("PA"),
            tags$td(167),
            tags$td(1990),
            tags$td(161)
          ),
          tags$tr(
            tags$td("PB1"),
            tags$td(113),
            tags$td(2165),
            tags$td(37)
          ),
          tags$tr(
            tags$td("PB2"),
            tags$td(109),
            tags$td(2152),
            tags$td(69)
          ),
          tags$tr(
            tags$td("PB2"),
            tags$td(163),
            tags$td(2152),
            tags$td(73)
          )
        )
      )

    ),
    box(
      title="What are 'direct repeats'?",
      width=12,
      "The sequence of the RNA before the starting point and end point of the",
      "deletion site can be the same (is repeated). This phenomena is",
      "described in literature a 'direct (sequence) repeat. Direct repeats",
      "can be of different length and are disscussed to be a driving factor",
      "in the generation of DelVGs.",
      br(),
      tags$img(src="direct_repeats.png"),
      br(),
      "The actual start and end point of a DelVG that has a direct repeat",
      "longer than 0 can not be determined certainly. In the given image a",
      "direct repeat of length 3 is depicted. There are four possible ways on",
      "how it could have been created during the experiment. If n is the",
      "length of the direct repeat, there are always n+1 options on how it",
      "could have been created."
    ),
    box(
      title="Contact information",
      width=12,
      "To get into contact with the developement team open a new",
      tags$a(href="https://github.com/LohmannJens/DIP-DSA/issues",
        target="_blank",
        "issue"
      ),
      "on GitHub. There you can get help with more detailed questions and",
      "come up with new ideas and features for the application."
    )
  )
)
