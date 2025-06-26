classdef globalsets
    properties (Constant)
        %% Cores per job
        NumberCoresStack = 18;
        NumberCoresCOLD  = 33;
        token = 'eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6ImxpbWVtZXJ0aHkiLCJleHAiOjE3MzkxNjMwMjQsImlhdCI6MTczMzk3OTAyNCwiaXNzIjoiaHR0cHM6Ly91cnMuZWFydGhkYXRhLm5hc2EuZ292IiwiaWRlbnRpdHlfcHJvdmlkZXIiOiJlZGxfb3BzIiwiYXNzdXJhbmNlX2xldmVsIjoyfQ.ELXYs5sEkqkzr3ib_E-IhfZp9BUhRNOFKVhNYJPWhL9KRIj8jXLoFr-FlOIIyFeRI5GqVNxCzkg4ihqo5ofmfYQKbJ7iwwAlKz8kspdzu1bhew8TMWovw8TQbjhasDlrmZCpTee9RcoII-RK2B4BQoRzbhd_bPBsyuarLikpGAR22tuqBMxU1_kpgMi_YAClvUfaQtS3yi7_j69s3Nf0g8ttqfYY6rVvGANilRKom14BKclc-veP9XZlW6tpRsl9MIh8sVr7nDRNLB1yH89LTHguv_uXL-K4hYQz6mMNGMnpCSUR3uCg8O1NnM6FLiZHOGdE2c1YWdL-l0pYSpC__g';
        
        %% Parameters
        Years = 2013:2023;
        DateStart = 2013001;
        DateEnd = 2023365;
        NTLCollectionVersion = '01';
        FolderpathStackCode = '/home/til19015/GlobalNTLAnalyze/StackData/';

        %% Paths % where '%s' will be replaced by NTL tile's name
        FolderpathNTLRaw = '/shared/cn449/VIIRS_NTL/BlackMarble_Raw/%s/';
        FolderpathNTLBRDF = '/shared/cn449/VIIRS_NTL/BlackMarble_BRDFcorrected/%s/';

        %% where '%s' will be replaced by NTL tile's name
        FolderpathStack = '/shared/zhulab/Tian/Stacked_new/%s/';

        %% Working Folders
        FolderCode = 'Code';
        FolderStack = 'stackdata';
    end
end

