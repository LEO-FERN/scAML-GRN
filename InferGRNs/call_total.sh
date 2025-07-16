echo "Network Inference..."

Rscript NetworkImputation&Inference.R --input /home/leandro/Data/Progenitor_Datasets --output /home/leandro/Data/Final_Progenitor_Net

echo "Progenitor Complete!"

Rscript NetworkImputation&Inference.R --input /home/leandro/Data/Monocyte_Datasets --output /home/leandro/Data/Final_Monocyte_Net

echo "Monocyte Complete!"

Rscript NetworkImputation&Inference.R --input /home/leandro/Data/Dendritic_Datasets --output /home/leandro/Data/Final_Dendritic_Net

echo "Dendritic Complete!"
echo "Networks Complete"
