function [res] = coreluacua (diskruptor1,diskruptor2)
res=norm(diskruptor1-diskruptor2)/(norm(diskruptor1)+1);

    