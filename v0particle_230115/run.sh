#o2-analysis-lf-lambdakzeroanalysis-mc  --configuration json://configuration.json -b   |
o2-analysis-mm-dndeta-hi --configuration json://configuration.json -b   |
o2-analysis-lf-lambdakzerobuilder - --configuration json://configuration.json  -b  |
o2-analysis-event-selection - --configuration json://configuration.json  -b  |
o2-analysis-multiplicity-table  --configuration json://configuration.json -b   --aod-file @input_data.txt |
o2-analysis-pid-tpc  --configuration json://configuration.json  -b  |
o2-analysis-timestamp  --configuration json://configuration.json  -b |
o2-analysis-track-propagation - --configuration json://configuration.json  -b |
o2-analysis-collision-converter - --configuration json://configuration.json  -b  |
o2-analysis-mm-track-propagation  -b --configuration json://configuration.json   --aod-file @input_data.txt |
o2-analysis-trackselection - --configuration json://configuration.json  -b  |
o2-analysis-centrality-table - --configuration json://configuration.json  -b   
