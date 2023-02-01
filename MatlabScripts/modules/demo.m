function demo()


OT1main('demo_keyfileOT1.txt','E:\Disser\Sets\Activity\sesq80.sdf');
OT2main('demo_keyfileOT2.txt','demo_distance.txt','demo_training_sample.sdf');
cluster_clusterdata('demo_keyfileOT1.txt', 'demo_y.txt')
cluster_clusterdata('demo_keyfileOT2.txt', 'demo_y.txt')
MGUAmd('demo_keyfileOT1.txt', 'demo_y.txt');
MGUAmd('demo_keyfileOT2.txt', 'demo_y.txt');
MGUAcluster('demo_keyfileOT2.txt');
generalCluster('demo_keyfileOT2.txt');
sampleDescriptionOT1('demo_test_keyfileOT1.txt', 'demo_keyfileOT1.txt', 'demo_test_sample.sdf');
sampleDescriptionOT2('demo_test_keyfileOT2.txt', 'demo_keyfileOT2.txt','demo_distance.txt', 'demo_test_sample.sdf');
predictionMD('demo_test_keyfileOT2.txt', 'demo_keyfileOT2.txt');
total_activityMD('demo_test_keyfileOT2.txt');
