create-batch --jobName checkabseta-re/deepIsoDY --fileList files_z.txt --cfg IsoConfig_cfg.py --maxFiles 10 --args 'Labels=1'
create-batch --jobName checkabseta-re/deepIsoQCD23 --fileList files_qcd-20to30.txt --cfg IsoConfig_cfg.py --maxFiles 10 --args 'Labels=0'
create-batch --jobName checkabseta-re/deepIsoQCD35 --fileList files_qcd-30to50.txt --cfg IsoConfig_cfg.py --maxFiles 10 --args 'Labels=0'
create-batch --jobName checkabseta-re/deepIsoQCD58 --fileList files_qcd-50to80.txt --cfg IsoConfig_cfg.py --maxFiles 10 --args 'Labels=0'
create-batch --jobName checkabseta-re/deepIsoQCD812 --fileList files_qcd-80to120.txt --cfg IsoConfig_cfg.py --maxFiles 10 --args 'Labels=0'
