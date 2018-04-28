%% find the actual numbers of occurrences for four tuples at interface with
%% two sites from each protein
function result = StatisticsWithCountNumbersForBigDatabase4nodes(protein1, EBindex1, resultpath)
%% statistics on datasets
num4nodetypes = 210;
fourcliqueType(:,:,1,1) =   [204    37    40    43    46    49    52
    37    38    99   102   105   108   111
    40    99    41   114   117   120   123
    43   102   114    44   126   129   132
    46   105   117   126    47   135   138
    49   108   120   129   135    50   141
    52   111   123   132   138   141    53];
fourcliqueType(:,:,2,1) =    [37    38    99   102   105   108   111
    38    36   100   103   106   109   112
    99   100   101     1     2     3     4
    102   103     1   104     5     6     7
    105   106     2     5   107     8     9
    108   109     3     6     8   110    10
    111   112     4     7     9    10   113];
fourcliqueType(:,:,3,1) =[ 40    99    41   114   117   120   123
    99   100   101     1     2     3     4
    41   101    39   115   118   121   124
    114     1   115   116    11    12    13
    117     2   118    11   119    14    15
    120     3   121    12    14   122    16
    123     4   124    13    15    16   125];
fourcliqueType(:,:,4,1) =[    43   102   114    44   126   129   132
    102   103     1   104     5     6     7
    114     1   115   116    11    12    13
    44   104   116    42   127   130   133
    126     5    11   127   128    17    18
    129     6    12   130    17   131    19
    132     7    13   133    18    19   134];
fourcliqueType(:,:,5,1) =[    46   105   117   126    47   135   138
    105   106     2     5   107     8     9
    117     2   118    11   119    14    15
    126     5    11   127   128    17    18
    47   107   119   128    45   136   139
    135     8    14    17   136   137    20
    138     9    15    18   139    20   140];
fourcliqueType(:,:,6,1) =[ 49   108   120   129   135    50   141
    108   109     3     6     8   110    10
    120     3   121    12    14   122    16
    129     6    12   130    17   131    19
    135     8    14    17   136   137    20
    50   110   122   131   137    48   142
    141    10    16    19    20   142   143];
fourcliqueType(:,:,7,1) =[    52   111   123   132   138   141    53
    111   112     4     7     9    10   113
    123     4   124    13    15    16   125
    132     7    13   133    18    19   134
    138     9    15    18   139    20   140
    141    10    16    19    20   142   143
    53   113   125   134   140   143    51];
fourcliqueType(:,:,1,2) =[    37    38    99   102   105   108   111
    38    36   100   103   106   109   112
    99   100   101     1     2     3     4
    102   103     1   104     5     6     7
    105   106     2     5   107     8     9
    108   109     3     6     8   110    10
    111   112     4     7     9    10   113];
fourcliqueType(:,:,2,2) =[    38    36   100   103   106   109   112
    36   205    55    58    61    64    67
    100    55    56   144   147   150   153
    103    58   144    59   156   159   162
    106    61   147   156    62   165   168
    109    64   150   159   165    65   171
    112    67   153   162   168   171    68];
fourcliqueType(:,:,3,2) =[    99   100   101     1     2     3     4
    100    55    56   144   147   150   153
    101    56    54   145   148   151   154
    1   144   145   146    21    22    23
    2   147   148    21   149    24    25
    3   150   151    22    24   152    26
    4   153   154    23    25    26   155];
fourcliqueType(:,:,4,2) =[   102   103     1   104     5     6     7
    103    58   144    59   156   159   162
    1   144   145   146    21    22    23
    104    59   146    57   157   160   163
    5   156    21   157   158    27    28
    6   159    22   160    27   161    29
    7   162    23   163    28    29   164];
fourcliqueType(:,:,5,2) =[   105   106     2     5   107     8     9
    106    61   147   156    62   165   168
    2   147   148    21   149    24    25
    5   156    21   157   158    27    28
    107    62   149   158    60   166   169
    8   165    24    27   166   167    30
    9   168    25    28   169    30   170];
fourcliqueType(:,:,6,2) =[   108   109     3     6     8   110    10
    109    64   150   159   165    65   171
    3   150   151    22    24   152    26
    6   159    22   160    27   161    29
    8   165    24    27   166   167    30
    110    65   152   161   167    63   172
    10   171    26    29    30   172   173];
fourcliqueType(:,:,7,2) =[   111   112     4     7     9    10   113
    112    67   153   162   168   171    68
    4   153   154    23    25    26   155
    7   162    23   163    28    29   164
    9   168    25    28   169    30   170
    10   171    26    29    30   172   173
    113    68   155   164   170   173    66];
fourcliqueType(:,:,1,3) =[    40    99    41   114   117   120   123
    99   100   101     1     2     3     4
    41   101    39   115   118   121   124
    114     1   115   116    11    12    13
    117     2   118    11   119    14    15
    120     3   121    12    14   122    16
    123     4   124    13    15    16   125];
fourcliqueType(:,:,2,3) =[    99   100   101     1     2     3     4
    100    55    56   144   147   150   153
    101    56    54   145   148   151   154
    1   144   145   146    21    22    23
    2   147   148    21   149    24    25
    3   150   151    22    24   152    26
    4   153   154    23    25    26   155];
fourcliqueType(:,:,3,3) =[    41   101    39   115   118   121   124
    101    56    54   145   148   151   154
    39    54   206    70    73    76    79
    115   145    70    71   174   177   180
    118   148    73   174    74   183   186
    121   151    76   177   183    77   189
    124   154    79   180   186   189    80];
fourcliqueType(:,:,4,3) =[   114     1   115   116    11    12    13
    1   144   145   146    21    22    23
    115   145    70    71   174   177   180
    116   146    71    69   175   178   181
    11    21   174   175   176    31    32
    12    22   177   178    31   179    33
    13    23   180   181    32    33   182];
fourcliqueType(:,:,5,3) =[   117     2   118    11   119    14    15
    2   147   148    21   149    24    25
    118   148    73   174    74   183   186
    11    21   174   175   176    31    32
    119   149    74   176    72   184   187
    14    24   183    31   184   185    34
    15    25   186    32   187    34   188];
fourcliqueType(:,:,6,3) =[   120     3   121    12    14   122    16
    3   150   151    22    24   152    26
    121   151    76   177   183    77   189
    12    22   177   178    31   179    33
    14    24   183    31   184   185    34
    122   152    77   179   185    75   190
    16    26   189    33    34   190   191];
fourcliqueType(:,:,7,3) =[   123     4   124    13    15    16   125
    4   153   154    23    25    26   155
    124   154    79   180   186   189    80
    13    23   180   181    32    33   182
    15    25   186    32   187    34   188
    16    26   189    33    34   190   191
    125   155    80   182   188   191    78];
fourcliqueType(:,:,1,4) =[    43   102   114    44   126   129   132
    102   103     1   104     5     6     7
    114     1   115   116    11    12    13
    44   104   116    42   127   130   133
    126     5    11   127   128    17    18
    129     6    12   130    17   131    19
    132     7    13   133    18    19   134];
fourcliqueType(:,:,2,4) =[   102   103     1   104     5     6     7
    103    58   144    59   156   159   162
    1   144   145   146    21    22    23
    104    59   146    57   157   160   163
    5   156    21   157   158    27    28
    6   159    22   160    27   161    29
    7   162    23   163    28    29   164];
fourcliqueType(:,:,3,4) =[   114     1   115   116    11    12    13
    1   144   145   146    21    22    23
    115   145    70    71   174   177   180
    116   146    71    69   175   178   181
    11    21   174   175   176    31    32
    12    22   177   178    31   179    33
    13    23   180   181    32    33   182];
fourcliqueType(:,:,4,4) =[    44   104   116    42   127   130   133
    104    59   146    57   157   160   163
    116   146    71    69   175   178   181
    42    57    69   207    82    85    88
    127   157   175    82    83   192   195
    130   160   178    85   192    86   198
    133   163   181    88   195   198    89];
fourcliqueType(:,:,5,4) =[   126     5    11   127   128    17    18
    5   156    21   157   158    27    28
    11    21   174   175   176    31    32
    127   157   175    82    83   192   195
    128   158   176    83    81   193   196
    17    27    31   192   193   194    35
    18    28    32   195   196    35   197];
fourcliqueType(:,:,6,4) =[   129     6    12   130    17   131    19
    6   159    22   160    27   161    29
    12    22   177   178    31   179    33
    130   160   178    85   192    86   198
    17    27    31   192   193   194    35
    131   161   179    86   194    84   199
    19    29    33   198    35   199   200];
fourcliqueType(:,:,7,4) =[   132     7    13   133    18    19   134
    7   162    23   163    28    29   164
    13    23   180   181    32    33   182
    133   163   181    88   195   198    89
    18    28    32   195   196    35   197
    19    29    33   198    35   199   200
    134   164   182    89   197   200    87];
fourcliqueType(:,:,1,5) =[    46   105   117   126    47   135   138
    105   106     2     5   107     8     9
    117     2   118    11   119    14    15
    126     5    11   127   128    17    18
    47   107   119   128    45   136   139
    135     8    14    17   136   137    20
    138     9    15    18   139    20   140];
fourcliqueType(:,:,2,5) =[   105   106     2     5   107     8     9
    106    61   147   156    62   165   168
    2   147   148    21   149    24    25
    5   156    21   157   158    27    28
    107    62   149   158    60   166   169
    8   165    24    27   166   167    30
    9   168    25    28   169    30   170];
fourcliqueType(:,:,3,5) =[   117     2   118    11   119    14    15
    2   147   148    21   149    24    25
    118   148    73   174    74   183   186
    11    21   174   175   176    31    32
    119   149    74   176    72   184   187
    14    24   183    31   184   185    34
    15    25   186    32   187    34   188];
fourcliqueType(:,:,4,5) =[   126     5    11   127   128    17    18
    5   156    21   157   158    27    28
    11    21   174   175   176    31    32
    127   157   175    82    83   192   195
    128   158   176    83    81   193   196
    17    27    31   192   193   194    35
    18    28    32   195   196    35   197];
fourcliqueType(:,:,5,5) =[    47   107   119   128    45   136   139
    107    62   149   158    60   166   169
    119   149    74   176    72   184   187
    128   158   176    83    81   193   196
    45    60    72    81   208    91    94
    136   166   184   193    91    92   201
    139   169   187   196    94   201    95];
fourcliqueType(:,:,6,5) =[   135     8    14    17   136   137    20
    8   165    24    27   166   167    30
    14    24   183    31   184   185    34
    17    27    31   192   193   194    35
    136   166   184   193    91    92   201
    137   167   185   194    92    90   202
    20    30    34    35   201   202   203];
fourcliqueType(:,:,7,5) =[   138     9    15    18   139    20   140
    9   168    25    28   169    30   170
    15    25   186    32   187    34   188
    18    28    32   195   196    35   197
    139   169   187   196    94   201    95
    20    30    34    35   201   202   203
    140   170   188   197    95   203    93];
fourcliqueType(:,:,1,6) =[    49   108   120   129   135    50   141
    108   109     3     6     8   110    10
    120     3   121    12    14   122    16
    129     6    12   130    17   131    19
    135     8    14    17   136   137    20
    50   110   122   131   137    48   142
    141    10    16    19    20   142   143];
fourcliqueType(:,:,2,6) =[   108   109     3     6     8   110    10
    109    64   150   159   165    65   171
    3   150   151    22    24   152    26
    6   159    22   160    27   161    29
    8   165    24    27   166   167    30
    110    65   152   161   167    63   172
    10   171    26    29    30   172   173];
fourcliqueType(:,:,3,6) =[   120     3   121    12    14   122    16
    3   150   151    22    24   152    26
    121   151    76   177   183    77   189
    12    22   177   178    31   179    33
    14    24   183    31   184   185    34
    122   152    77   179   185    75   190
    16    26   189    33    34   190   191];
fourcliqueType(:,:,4,6) =[   129     6    12   130    17   131    19
    6   159    22   160    27   161    29
    12    22   177   178    31   179    33
    130   160   178    85   192    86   198
    17    27    31   192   193   194    35
    131   161   179    86   194    84   199
    19    29    33   198    35   199   200];
fourcliqueType(:,:,5,6) =[   135     8    14    17   136   137    20
    8   165    24    27   166   167    30
    14    24   183    31   184   185    34
    17    27    31   192   193   194    35
    136   166   184   193    91    92   201
    137   167   185   194    92    90   202
    20    30    34    35   201   202   203];
fourcliqueType(:,:,6,6) =[    50   110   122   131   137    48   142
    110    65   152   161   167    63   172
    122   152    77   179   185    75   190
    131   161   179    86   194    84   199
    137   167   185   194    92    90   202
    48    63    75    84    90   209    97
    142   172   190   199   202    97    98];
fourcliqueType(:,:,7,6) =[   141    10    16    19    20   142   143
    10   171    26    29    30   172   173
    16    26   189    33    34   190   191
    19    29    33   198    35   199   200
    20    30    34    35   201   202   203
    142   172   190   199   202    97    98
    143   173   191   200   203    98    96];
fourcliqueType(:,:,1,7) =[    52   111   123   132   138   141    53
    111   112     4     7     9    10   113
    123     4   124    13    15    16   125
    132     7    13   133    18    19   134
    138     9    15    18   139    20   140
    141    10    16    19    20   142   143
    53   113   125   134   140   143    51];
fourcliqueType(:,:,2,7) =[   111   112     4     7     9    10   113
    112    67   153   162   168   171    68
    4   153   154    23    25    26   155
    7   162    23   163    28    29   164
    9   168    25    28   169    30   170
    10   171    26    29    30   172   173
    113    68   155   164   170   173    66];
fourcliqueType(:,:,3,7) =[   123     4   124    13    15    16   125
    4   153   154    23    25    26   155
    124   154    79   180   186   189    80
    13    23   180   181    32    33   182
    15    25   186    32   187    34   188
    16    26   189    33    34   190   191
    125   155    80   182   188   191    78];
fourcliqueType(:,:,4,7) =[   132     7    13   133    18    19   134
    7   162    23   163    28    29   164
    13    23   180   181    32    33   182
    133   163   181    88   195   198    89
    18    28    32   195   196    35   197
    19    29    33   198    35   199   200
    134   164   182    89   197   200    87];
fourcliqueType(:,:,5,7) =[   138     9    15    18   139    20   140
    9   168    25    28   169    30   170
    15    25   186    32   187    34   188
    18    28    32   195   196    35   197
    139   169   187   196    94   201    95
    20    30    34    35   201   202   203
    140   170   188   197    95   203    93];
fourcliqueType(:,:,6,7) =[   141    10    16    19    20   142   143
    10   171    26    29    30   172   173
    16    26   189    33    34   190   191
    19    29    33   198    35   199   200
    20    30    34    35   201   202   203
    142   172   190   199   202    97    98
    143   173   191   200   203    98    96];
fourcliqueType(:,:,7,7) =[    53   113   125   134   140   143    51
    113    68   155   164   170   173    66
    125   155    80   182   188   191    78
    134   164   182    89   197   200    87
    140   170   188   197    95   203    93
    143   173   191   200   203    98    96
    51    66    78    87    93    96   210];


if ~iscell(protein1)
    protein{1} = protein1;   
else
    protein = protein1;
end
if ~iscell(EBindex1)
    EBindex{1} = EBindex1;
else
    EBindex = EBindex1;
end

aminonames = {'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'A'};

%% find all the interdomain contact sites and interdomain contact pairs for
% % protein.aln.nogap
% 'find contactSites and contactPairs ...'
clear contactSite; clear contactPair; clear noncontactSite;
for pri = 1 : size(protein,2)
    fiddata =  fopen([resultpath, 'contactInformation', protein{pri}.names, '.mat']);
    if fiddata == -1        
        contactPair{pri} = [];
        clear location;
        location = get_locationMSA_from_resnum(protein{pri}.aln.nogap.locationMSA,...
            protein{pri}.aln.nogap.residue, protein{pri}.aln.contact_pairs_residue.number,...
            protein{pri}.aln.contact_pairs_residue.chain);
        % record inter-domain contact pairs
        for i = 1 : size(location,1)
            if location(i,1) * location(i,2)  ~= 0
                if EBindex{pri}(location(i,1)) == 2 && EBindex{pri}(location(i,2)) == 2
                    % interdomain contact pairs that are actually on the
                    % surface
                    if location(i,1) <= protein{pri}.aln.nogap.seprate_site && location(i,2) > protein{pri}.aln.nogap.seprate_site ...
                            || location(i,2) <= protein{pri}.aln.nogap.seprate_site && location(i,1) > protein{pri}.aln.nogap.seprate_site
                        contactPair{pri} = [contactPair{pri}; location(i,:)];
                    end
                end
            end
        end
        % get the contact sites from the inter-domain contact pairs on the
        % surface
        contactSite{pri} = [];
        noncontactSite{pri} = [];
        for i = 1 : size(contactPair{pri},1)
            for j = 1 : 2
                if ~isin(contactPair{pri}(i,j), contactSite{pri})
                    contactSite{pri} = [contactSite{pri}, contactPair{pri}(i,j)];
                end
            end
        end
        for i = 1 : protein{pri}.aln.nogap.length
            if ~isin(i,contactSite{pri})
                noncontactSite{pri} = [noncontactSite{pri}, i];
            end
        end  
        clear temp; 
        temp.contactSite = contactSite{pri};
        temp.noncontactSite = noncontactSite{pri};
        temp.contactPair = contactPair{pri};
        save([resultpath, 'contactInformation', protein{pri}.names, '.mat'], 'temp');
        clear temp;
    else
        fclose(fiddata);
        clear temp;
        temp = load([resultpath, 'contactInformation', protein{pri}.names, '.mat']);
        contactSite{pri} = temp.temp.contactSite;
        noncontactSite{pri} = temp.temp.noncontactSite;
        contactPair{pri} = temp.temp.contactPair;
        clear temp;
    end
end

%% map amino acids into indexes
% 'map amino acids into indexes ...'
clear IndexSeq;
for pri = 1 : size(protein,2)  
    % for training data
    fiddata =  fopen([resultpath, 'mappedindex', protein{pri}.names, '.mat']);
    if fiddata == -1
        clear mappedindex;
        indexPDB = 1;
        for i = 1 : protein{pri}.aln.nogap.length
            mappedindex(i) = MapResiduesToIndex(protein{pri}.aln.nogap.seq(indexPDB).Sequence(i), aminonames);
        end
        IndexSeq{pri} = mappedindex;
        save([resultpath, 'mappedindex', protein{pri}.names, '.mat'], 'mappedindex');
    else
        fclose(fiddata);
        clear temp;
        temp = load([resultpath, 'mappedindex', protein{pri}.names, '.mat']);
        IndexSeq{pri} = temp.mappedindex;
        clear temp;
    end
end

%%  mapinto20AA
clear IndexSeq20;
for pri = 1 : size(protein,2)
    fiddata = fopen([resultpath, 'IndexMSA20AA', protein{pri}.names,'.mat']);
    if fiddata == -1
        for i = 1 : protein{pri}.aln.nogap.length
            mapped_index(i) = MapResiduesToIndex20(protein{pri}.aln.nogap.seq(1).Sequence(i), aminonames);
        end
        IndexSeq20{pri} = mapped_index;
        save([resultpath, 'IndexMSA20AA', protein{pri}.names,'.mat'], 'mapped_index');
    else
        fclose(fiddata);
        clear tempdata;
        tempdata = load([resultpath, 'IndexMSA20AA', protein{pri}.names,'.mat']);
        IndexSeq20{pri} = tempdata.mapped_index;
        clear tempdata;
    end
end


%% structure graph
% 'form structure graph ...'
for pri = 1 : size(protein,2)  
    % for training data
    fiddata = fopen([resultpath, 'structuregraph', protein{pri}.names, '.mat']);
    if fiddata == -1
        clear graph;
        graph = make_structure_graph(protein{pri});
        structuregraph{pri} = graph;        
        save([resultpath, 'structuregraph', protein{pri}.names, '.mat'], 'graph');
    else
        fclose(fiddata);
        clear temp;
        temp = load([resultpath, 'structuregraph', protein{pri}.names, '.mat']);
        structuregraph{pri} = temp.graph;
        clear temp;
    end
end

%% interdomain contact :
count_con_clique = zeros(num4nodetypes, 11);   %
for pri = 1 : size(protein,2)
    % separate_site
    seprate_site = protein{pri}.aln.nogap.seprate_site;
    % contact sites on two domains
    sitesondomain{1} = [];
    sitesondomain{2} = [];
    clear structuregraph1;
    structuregraph1 = structuregraph{pri};
    for i = 1 : size(contactSite{pri}, 2)
        if contactSite{pri}(i) <= seprate_site
            sitesondomain{1} = [sitesondomain{1}, contactSite{pri}(i)];
        else
            sitesondomain{2} = [sitesondomain{2}, contactSite{pri}(i)];
        end
    end
    sitesondomain1 = zeros(1, 2000);      s1 = size(sitesondomain{1}, 2);  sitesondomain1(1, 1: s1) = sitesondomain{1};
    sitesondomain2 = zeros(1,2000);       s2 = size(sitesondomain{2}, 2);  sitesondomain2(1, 1: s2) = sitesondomain{2};
    indexseqall = zeros(1,4000);              indexseqall(1, 1:protein{pri}.aln.nogap.length) = IndexSeq{pri} ;
    structuregraph2 = zeros(2000,2000);  structuregraph2(1:protein{pri}.aln.nogap.length, 1:protein{pri}.aln.nogap.length) = structuregraph1;
    count_con_clique = CliqueCounting4forcon(sitesondomain1, s1,  sitesondomain2, s2, count_con_clique, fourcliqueType, indexseqall, structuregraph2);    
end
%% commented on 12th May
% %% interdomain surface background:
% count_non_clique = zeros(num4nodetypes, 6);   %
% for pri = 1 : size(protein,2)
%     % separate_site
%     seprate_site = protein{pri}.aln.nogap.seprate_site;
%     % contact sites on two domains
%     sitesondomain{1} = [];
%     sitesondomain{2} = [];
%     % surface sites only
%     location = find(EBindex{pri}(1,:) == 2);     
%     % only remain those non-contact sites
%     for i = 1 : size(location,2)
%         if  location(i) <= seprate_site
%             sitesondomain{1} = [sitesondomain{1}, location(i)];
%         else
%             sitesondomain{2} = [sitesondomain{2}, location(i)];
%         end
%     end    
%     sitesondomain1 = zeros(1, 2000);      s1 = size(sitesondomain{1}, 2);  sitesondomain1(1, 1: s1) = sitesondomain{1};
%     sitesondomain2 = zeros(1,2000);       s2 = size(sitesondomain{2}, 2);  sitesondomain2(1, 1: s2) = sitesondomain{2};
%     indexseqall = zeros(1,4000);              indexseqall(1, 1:protein{pri}.aln.nogap.length) = IndexSeq{pri} ;
%     count_non_clique = CliqueCounting4forall(sitesondomain1, s1, sitesondomain2, s2, count_non_clique, fourcliqueType, indexseqall);
% end
%% end of commented on 12th May

%    clear structuregraph1;
%    structuregraph1 = structuregraph{pri};
%     % exclude the contact sites
%     location = find(EBindex{pri}(1,:) == 2);     
%     for i = 1 : size(contactSite{pri},2)
%         location(location==contactSite{pri}(i)) = [];
% %     end
%     % only remain those non-contact sites
%     for i = 1 : size(location,2)
%         if  location(i) <= seprate_site
%             sitesondomain{1} = [sitesondomain{1}, location(i)];
%         else
%             sitesondomain{2} = [sitesondomain{2}, location(i)];
%         end
%     end  
% % pick up two sites from each domain to form the 4-clique candidate
%     for i = 1 : size(sitesondomain{1},2)-1
%         for j = i+1 : size(sitesondomain{1},2)
%             for p = 1 : size(sitesondomain{2},2)-1
%                 for q = p+1 : size(sitesondomain{2},2)
%                     cliquetyep = ClassifyTypes(IndexSeq{pri}([sitesondomain{1}(i),...
%                         sitesondomain{1}(j), sitesondomain{2}(p), sitesondomain{2}(q)]));
%                     if cliquetyep > 0
%                         count_all_clique(cliquetyep,1) = count_all_clique(cliquetyep,1) + 1;
%                     end
%                 end
%             end
%         end
%     end

% result.count_non_clique = count_non_clique;
result.count_con_clique = count_con_clique;
result.contactSite = contactSite;
result.contactPair = contactPair;


%% get residue distances if they have structure info. assign distance to
% nogap residues
% protein.nogap
% since on the big database there is no protein.aln.nogap.structure 
% every residue has structure info, so the location1 is actually the
% locatioin in the distance matrix
function distance = get_distance(location1, location2, protein)
% residue
residue1 = get_resnum_from_locationMSA(protein.aln.nogap.locationMSA, protein.aln.nogap.residue, location1);
residue2 = get_resnum_from_locationMSA(protein.aln.nogap.locationMSA, protein.aln.nogap.residue, location2);
location_structure1 = get_locationMSA_from_resnum(protein.aln.nogap.locationMSA, protein.aln.nogap.residue, ...
    residue1.number, residue1.chain);
location_structure2 = get_locationMSA_from_resnum(protein.aln.nogap.locationMSA, protein.aln.nogap.residue, ...
    residue2.number, residue2.chain);
if location_structure1 * location_structure2 ~= 0
    distance = protein.aln.nogap.distances(location_structure1, location_structure2);
else
    distance = 100000;
end


%% make distances graph for protein, if atom-atom distances smaller than
% 4.5, a contact is identified
function graph = make_structure_graph(protein)
graph = zeros(protein.aln.nogap.length);
for i = 1 : protein.aln.nogap.length-1
    for j = i+1 : protein.aln.nogap.length
        if get_distance(i,j,protein) < 4.5
            graph(i,j) = 1;
            graph(j,i) = graph(i,j);
        end
    end
end

%% find the two sites from the same domain and the sites in the other domain
function [a_intra,b_intra,c_inter] = find_intra_pair(i,j,k, seprate_site)
domain1 = []; domain2 = [];
if i <= seprate_site
    domain1 = [domain1, i];
else
    domain2 = [domain2, i];
end
if j <= seprate_site
    domain1 = [domain1, j];
else
    domain2 = [domain2, j];
end
if k <= seprate_site
    domain1 = [domain1, k];
else
    domain2 = [domain2, k];
end
if size(domain1,2) > 1
    a_intra = domain1(1);
    b_intra = domain1(2);
    if size(domain2,2) > 1
        'error in triangle define'
    else
        c_inter = domain2(1);
    end
else
    if size(domain2,2) > 1
        a_intra = domain2(1);
        b_intra = domain2(2);
    else
        'error in triangle define'
    end
    c_inter = domain1(1);
end

    
% % is there any node with degree of three
% function bool = getdegreeofnode(structuregraph,sitesondomain, i,j,p,q)
% 
% summ1 =  structuregraph(sitesondomain{1}(i), sitesondomain{1}(j)) ...
%     + structuregraph(sitesondomain{1}(i), sitesondomain{2}(p)) ...
%     + structuregraph(sitesondomain{1}(i), sitesondomain{2}(q));
% 
% summ2 =  structuregraph(sitesondomain{1}(j), sitesondomain{1}(i)) ...
%     + structuregraph(sitesondomain{1}(j), sitesondomain{2}(p)) ...
%     + structuregraph(sitesondomain{1}(j), sitesondomain{2}(q));
% 
% summ3 =  structuregraph(sitesondomain{2}(p), sitesondomain{1}(i)) ...
%     + structuregraph(sitesondomain{2}(p), sitesondomain{1}(j)) ...
%     + structuregraph(sitesondomain{2}(p), sitesondomain{2}(q));
% 
% summ4 =  structuregraph(sitesondomain{2}(q), sitesondomain{1}(i)) ...
%     + structuregraph(sitesondomain{2}(q), sitesondomain{1}(j)) ...
%     + structuregraph(sitesondomain{2}(q), sitesondomain{2}(p));
% 
% if summ1 == 3 || summ2 == 3 || summ3 == 3 || summ4 == 3
%     bool = 1;
% else
%     bool = 0;
% end





