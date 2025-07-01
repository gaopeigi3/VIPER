
java -Xmx5G -jar aracne.jar -e ../venetoclax_expression.tsv -o ../output --tfs ../tfs.txt --pvalue 1E-8 --seed 1 --calculateThreshold

for i in {1..100}
do
  java -Xmx5G -jar aracne.jar -e ../venetoclax_expression.tsv -o ../output --tfs ../tfs.txt --pvalue 1E-8 --seed $i
done

java -Xmx5G -jar aracne.jar -o ../output --consolidate
