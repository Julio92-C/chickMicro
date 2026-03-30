import PyPDF2


merger = PyPDF2.PdfFileMerger()

file_names = ["./kraken2_taxaCount.pdf",
              "./kraken2_vennDiagram.pdf",
              "./kraken2_relativeAbundance.pdf",
	      "./abriKraken2_taxaCount.pdf",
	      "./abrikraken2_relativeAbundance_legend.pdf",
              "./abriKraken2_ARGsCount.pdf",
	      "./abriKraken2_MGEsCount.pdf",
              "./GEsProfile_pieChart.pdf",
              "./ARGsProfile_chordDiagram.pdf",
              "./GEs_vennDiagram.pdf",
              "./ARGsCount_violinPlot.pdf",
              "./VFsCount_violinPlot.pdf",
	      "./VFsProfile_SankeyNetwork.pdf,
              "./GEsProfile_heatMap.pdf",
              ]


for file_name in file_names:
    merger.append(file_name)

merger.write("./HosMicro_Summary.pdf")
merger.close()
