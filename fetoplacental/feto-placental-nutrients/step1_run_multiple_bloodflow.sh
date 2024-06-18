list=(20.0 40.0 50.0 60.0 80.0 100.0 150.0 200.0 250.0 300.0 350.0 400.0 450.0 500.0)
for i in ${list[@]}
do
   echo $i
   python run-bloodflow.py $i 0.0  > screen_$i.txt
   mv screen_$i.txt output_flowbc_${i}ml_min_0.0
   #tar -cvf output_${i}ml_min.tar  output_flowbc_${i}ml_min_0.0
   #gzip output_${i}ml_min.tar
   #rm -rf output_flowbc_${i}ml_min_0.0
done

mkdir output-normal
mv output_flowbc_* output-normal

