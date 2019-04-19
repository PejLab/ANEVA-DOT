#' Vg_GTEx_v7
#'
#' A dataset containing the estimated variance in gene expression obtained
#' from the GTEx project. Use the command View(VG_GTEx_v7) to have a look at
#' the data.
#'
#' The Genotype-Tissue Expression (GTEx) Project was supported by the Common
#' Fund of the Office of the Director of the National Institutes of Health,
#' and by NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. The data used for the
#' analyses described in this manuscript were obtained from the GTEx Portal
#' on 02/17/19.
#'
#' @format A data frame with 14091 rows and 51 variables:
#' \describe{
#'   \item{IDs}{Ensembl ID}
#'   \item{BRNACC}{}
#'   \item{ARTCRN}{}
#'   \item{BRNHPT}{}
#'   \item{HRTLV}{}
#'   \item{ADPVSC}{}
#'   \item{CLNTRN}{}
#'   \item{PNCREAS}{}
#'   \item{TESTIS}{}
#'   \item{LUNG}{}
#'   \item{ADRNLG}{}
#'   \item{BREAST}{}
#'   \item{BLDDER}{}
#'   \item{SKINS}{}
#'   \item{SPLEEN}{}
#'   \item{BRNAMY}{}
#'   \item{BRNCHB}{}
#'   \item{BRNCDT}{}
#'   \item{BRNCTXB}{}
#'   \item{KDNCTX}{}
#'   \item{PRSTTE}{}
#'   \item{BRNHPP}{}
#'   \item{UTERUS}{}
#'   \item{CLNSGM}{}
#'   \item{ESPMCS}{}
#'   \item{PTTARY}{}
#'   \item{SKINNS}{}
#'   \item{ESPGEJ}{}
#'   \item{ARTTBL}{}
#'   \item{LCL}{}
#'   \item{MSCLSK}{}
#'   \item{SNTTRM}{}
#'   \item{BRNSNG}{}
#'   \item{THYROID}{}
#'   \item{STMACH}{}
#'   \item{FIBRBLS}{}
#'   \item{VAGINA}{}
#'   \item{BRNPTM}{}
#'   \item{LIVER}{}
#'   \item{BRNNCC}{}
#'   \item{BRNSPC}{}
#'   \item{WHLBLD}{}
#'   \item{ADPSBQ}{}
#'   \item{OVARY}{}
#'   \item{SLVRYG}{}
#'   \item{HRTAA}{}
#'   \item{BRNCHA}{}
#'   \item{ESPMSL}{}
#'   \item{BRNCTXA}{}
#'   \item{ARTAORT}{}
#'   \item{NERVET}{}
#' }
#' @source \url{https://gtexportal.org/home/}
"Vg_GTEx_v7"

#'Sample ASE data
#'
#'A sample of Allele Specific Expression (ASE) data from a single tissue sample from a GTEx v7 donor
#'for use with the short tutorial on the package README.
#'The data was produced using phASEr \url{https://github.com/secastel/phaser}
#'
#'@format A data frame with 5460 rows and 6 variables:
#'\describe{
#'  \item{GENE_ID}{The Ensembl Gene ID}
#'  \item{TISSUE_ID}{The tissue ID from GTEx. (Skeletal Muscle in this case)}
#'  \item{REF_COUNT}{The number of RNA-seq reads from the highest expressed ASE-snp mapping to the reference allele of the gene}
#'  \item{ALT_COUNT}{The number of RNA-seq reads from the highest expressed ASE-snp mapping to the alternative allele of the gene}
#'  \item{TOTAL_COUNT}{The sum of REF_COUNT and ALT_COUNT}
#'  \item{NULL_RATIO}{From phASEr, the expected value of ALT_COUNT/(ALT_COUNT + REF_COUNT)}
#'}
"sample_ASE"
