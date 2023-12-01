# grn
sudo docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like:/home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like     aertslab/pyscenic:0.12.1 pyscenic grn         --num_workers 25         -o /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/adj.retinoma_like.tsv         /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/matrix_retinoma_like.csv         /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/allTFs_hg38.txt
sudo docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased:/home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased     aertslab/pyscenic:0.12.1 pyscenic grn         --num_workers 25         -o /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/adj.MKI67_photoreceptorness_decreased.tsv         /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/matrix_MKI67_photoreceptorness_decreased.csv         /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/allTFs_hg38.txt
sudo docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like:/home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like     aertslab/pyscenic:0.12.1 pyscenic grn         --num_workers 25         -o /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/adj.Cone_precursor_like.tsv         /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/matrix_Cone_precursor_like.csv         /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/allTFs_hg38.txt

# ctx
docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like:/home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like     aertslab/pyscenic:0.12.1 pyscenic ctx         /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/adj.retinoma_like.tsv         /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather         /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather         --annotations_fname /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl         --expression_mtx_fname /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/matrix_retinoma_like.csv         --mode custom_multiprocessing         --output /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/reg.retinoma_like.csv         --num_workers 25
docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased:/home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased     aertslab/pyscenic:0.12.1 pyscenic ctx         /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/adj.MKI67_photoreceptorness_decreased.tsv         /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather         /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather         --annotations_fname /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl         --expression_mtx_fname /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/matrix_MKI67_photoreceptorness_decreased.csv         --mode custom_multiprocessing         --output /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/reg.MKI67_photoreceptorness_decreased.csv         --num_workers 25
docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like:/home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like     aertslab/pyscenic:0.12.1 pyscenic ctx         /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/adj.Cone_precursor_like.tsv         /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather         /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather         --annotations_fname /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl         --expression_mtx_fname /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/matrix_Cone_precursor_like.csv         --mode custom_multiprocessing         --output /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/reg.Cone_precursor_like.csv         --num_workers 25

# aucell
docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like:/home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like     aertslab/pyscenic:0.12.1 pyscenic aucell         /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/matrix_retinoma_like.csv         /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/reg.retinoma_like.csv         -o /home/weihu/rb_revision_2_project/results/5_SCENIC/retinoma_like/auc_mtx.retinoma_like.csv         --num_workers 15
docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased:/home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased     aertslab/pyscenic:0.12.1 pyscenic aucell         /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/matrix_MKI67_photoreceptorness_decreased.csv         /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/reg.MKI67_photoreceptorness_decreased.csv         -o /home/weihu/rb_revision_2_project/results/5_SCENIC/MKI67_photoreceptorness_decreased/auc_mtx.MKI67_photoreceptorness_decreased.csv         --num_workers 15
docker run -d --rm     -v /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like:/home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like     aertslab/pyscenic:0.12.1 pyscenic aucell         /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/matrix_Cone_precursor_like.csv         /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/reg.Cone_precursor_like.csv         -o /home/weihu/rb_revision_2_project/results/5_SCENIC/Cone_precursor_like/auc_mtx.Cone_precursor_like.csv         --num_workers 15