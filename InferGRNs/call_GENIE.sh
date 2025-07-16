echo "GENIE3..."

Rscript GENIE3_Inference.R --input /home/leandro/Data/Final_Progenitor_Net --output /home/leandro/Data/Final_Progenitor_Net

echo "Progenitor Complete!"

Rscript GENIE3_Inference.R --input /home/leandro/Data/Final_Monocyte_Net --output /home/leandro/Data/Final_Monocyte_Net

echo "Monocyte Complete!"

Rscript GENIE3_Inference.R --input /home/leandro/Data/Final_Dendritic_Net --output /home/leandro/Data/Final_Dendritic_Net

echo "Dendritic Complete!"
