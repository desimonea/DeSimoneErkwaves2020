for %%f in (G:\H2A_ERKKTR_test\TGMM_hypo_eq_ch2\TGMMconfig\*.txt) do (

       cd nucleiChSvWshedPBC\Release
       start /wait ProcessStackBatchMulticore.exe %%f 0 1
       cd ../..
       cd Release
       start /wait TGMM.exe %%f 0 1
       cd ..
       del G:\H2A_ERKKTR_test\TGMM_hypo_eq_ch2\TGMM_hypo_eq_ch2\data\*.bin
        
)

