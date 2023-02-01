function all_keys(f,KeyFile, param1,param2)

warning off
markers={'___','__r','_b_','d__','_br','d_r','db_','dbr'};
length={'2','3','4'};
KeyStruct=readKeys(KeyFile);
for I=1:size(markers,2)
    markers{I}
    for J=1:size(length,2)
        length{J}
        KeyStruct.chain_length=length{J};
        KeyStruct.markers=markers{I};
        setKeys(KeyFile, KeyStruct);
        f(KeyFile, param1);
    end
end

warning on

