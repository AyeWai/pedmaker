
#One important thing to note here, is you need to add the @export tag above your function
#to indicate this function to be "exposed" to users to use.

#The #' @export syntax is actually an Roxygen tag . By doing this, this ensures that the load_mat() function gets added to the NAMESPACE (when you run devtools::document()) to indicate that it needs to be exposed.

#' Life Hist Data creator
#'
#' Creates a LifeHistData taking in parameters the link file between genetic individuals names
#' and caught individuals names and the differents way to calculate the individuals age
#' and transforms the F/M sex in 1/2 sex .
#'
#' @param genecap_file The genetic-catch names links file
#' @param low Age is calculated with the a lower threshold with low = TRUE
#' @param high Age is calculated with the a upper threshold with high = TRUE
#' @param moyenne The mean age is calculated whith moyenne = TRUE
#' @param age.supp The age attributed to the 6 or more years old individuals
#' @return A LifeHistData file to be used to by PediCre
#' @export
LHD <- function(genecap_file,low = FALSE,high = FALSE,moyenne=FALSE,age.supp = 8){
  sink('Log_LHD.txt', append = TRUE, split = TRUE )
  print(Sys.time())
  print(c(genecap_file,low,high,moyenne,age.supp))
  data=read_excel(genecap_file)
  data1=subset(data,(ani_cor %nin% c("494"))) ###494 répété !
  data=droplevels(data1)
  data$ani_cor=as.character(data$ani_cor)
  data=data[,'ani_cor']
  assign('Corresp_NoDouble',data, envir = globalenv())
  drv <- dbDriver("PostgreSQL")
  ## FULL join permet de garder toutes les lignes du tableau
  con <- dbConnect(drv, dbname=dbn, host=hst, user=usr, password=pw) ## Ã  complÃ©ter

  mydata<-dbGetQuery(con,"SELECT public.t_animal_ani.*
                   from public.t_animal_ani;")

  mydata2<-dbGetQuery(con,"SELECT *
                      from t_capture_cap where cap_id in (
                      select toto.cap_id from (SELECT min(cap_id) as cap_id, cap_ani_id, min(cap_date)
                      FROM public.t_capture_cap, t_animal_ani where cap_ani_id = ani_id group by cap_ani_id) as toto) order by cap_id")

  mydata=merge(mydata,mydata2,by.x="ani_id",by.y="cap_ani_id",all.x=T,all.y=T)

  mydata3<-dbGetQuery(con,"SELECT public.t_animal_ani.*,public.t_capture_cap.cap_id,public.t_capture_cap.cap_faon,
                      public.t_capture_cap.cap_bague,public.t_capture_cap.cap_age_corrige,public.t_capture_cap.cap_age,public.t_capture_cap.cap_annee_suivi
                      from public.t_animal_ani,public.t_capture_cap where cap_ani_id=ani_id;")

  y = names(mydata3)
  mydata = mydata[,y]

  assign('DB_data',mydata, envir = globalenv())
  dbDisconnect(con)
  #### On bosse toujours sur ani_etiq pour le nom de l'individus
  dv1<-subset(mydata,ani_etiq %in% data$ani_cor)

  #### Creation LifeHistData brut

  dv2 = data.frame(dv1$ani_etiq,dv1$ani_sexe,dv1$cap_age, dv1$cap_annee_suivi)
  assign('Raw_LifeHistData',dv2, envir = globalenv())

  ##Calcul colonne age estimé individus lors de la capture
  nom = ''
  for (x in 1:nrow(dv2)) {
    if(grepl("[0-9]-[0-9]",dv2$dv1.cap_age[x]) == T) {
      minimum = as.integer(substr(dv2$dv1.cap_age[x],1,1))
      maximum = as.integer(substr(dv2$dv1.cap_age[x],3,3))
      if (moyenne){
      dv2$mean_age[x] =as.integer(mean(minimum,maximum))+1
      nom = 'Mean'
      }
      else if(low){
        dv2$mean_age[x] =minimum+1
        nom = 'Low'
      }
      else{
        dv2$mean_age[x] =maximum+1
        nom = 'Max'
      }
      }
    else if(grepl("^[1-9]",dv2$dv1.cap_age[x]) == T){

      dv2$mean_age[x] = as.integer(substr(dv2$dv1.cap_age[x],1,1)) +1 #car dv2$dv1.cap_age contient un facteur, et directement nous renvoie si on appelle dv2$dv1.cap_age[x], le niveau de dv2$dv1.cap_age[x]et non sa valeur.
    }
    else {
      if (grepl("^[>=]+6",dv2$dv1.cap_age[x]) == T){
        dv2$mean_age[x] = age.supp
      }
      else if (grepl("^<",dv2$dv1.cap_age[x]) == T){
        dv2$mean_age[x] = 1
      }
    }
  }
  ## Calcul date naissance estimee

  dv2$ani_by = dv2$dv1.cap_annee_suivi - as.integer(dv2$mean_age)

  ##Transformation du sexe individus F/M -> 1/2

  for(x in 1:nrow(dv2)) {
    if(dv2$dv1.ani_sexe[x] == "F"){
      dv2$ani_sexe2[x] = 1
    }
    else {
      dv2$ani_sexe2[x] = 2
    }
  }

  ##Rename des colonnes pour plus de clarté
  names(dv2) = c('ani_etiq','ani_sexeFM','cap_age','cap_annee_suivi','ani_mean_age','ani_by','ani_sexe')

  ##Tableau données/historique vie final
  LifeHistData=dv2[,c("ani_etiq","ani_sexe","ani_by")]
  assign('LifeHistData',LifeHistData, envir = globalenv())
  dv3 = data.frame(dv1$ani_etiq,dv1$ani_sexe,dv2$ani_sexe,dv1$cap_age, dv2$ani_mean_age,dv1$cap_annee_suivi, dv2$ani_by)
  assign('Intermediary_LifeHistData',dv3,envir=globalenv())
  # assign('LifeHisData', LifeHisData, envir = globalenv())
  write.table(LifeHistData, file =  paste('LifeHistData',nom,".txt",sep=""), sep = '\t', row.names = FALSE, quote =FALSE)
  sink()
}

#' Pedigree Creator
#'
#' Creates a Pedigree taking in paramater a raw file containing a 0/1/2 genome, the link file between genetic individuals names
#' and caught individuals names.
#' It also take the Maximum Siblings number Iterations number to relatives assignment as MS1 to the Parent-Offspring assignment
#' and MS2 to the others relatives assignement.
#'
#' @param raw_file The raw file name
#' @param genecap_file The genetic-catch names links file
#' @param LHData LifeHistData file
#' @param MS1 First run maximum number of sibship clustering iterations
#' @return The population's SNP pedigree
#' @export
#'
PediCre <- function(raw_file,genecap_file,LHData,MS){
  bbtime = as.character(Sys.time())
  bbsep = ('\n---------------------------------------------------------------------------------------------------\n')
  bb = file('Log_PediCre.txt', open = 'at')
  cat(bbsep, file = bb)
  sink(bb, append = TRUE, type = 'message')
  # sink('logfile_Pedicre.txt', type = 'message',append = TRUE, split = TRUE )
  correspondance_genetcapture <- read_excel(genecap_file)
  test_sequoia = read.table(raw_file, sep = ' ', header = T)
  assign('RawSequencing_PlinkFile',test_sequoia, envir = globalenv())
  preGeno = merge(test_sequoia, correspondance_genetcapture, by.x = 'IID', by.y = 'ani_genet' )
  preGeno$IID = preGeno$ani_cor
  preGeno$FID = preGeno$ani_cor
  preGeno = preGeno[,1:406]

  write.table(preGeno, file = paste('A-',raw_file,sep = ''), sep = ' ', row.names = FALSE, quote =FALSE)

  Geno = GenoConvert(InFile = paste('A-',raw_file,sep = ''))
  assign('GenoM',Geno, envir = globalenv())
  XX = read.table(LHData,header = T, sep = '\t')
  assign('LifeHistData',XX, envir = globalenv())

  ParOUT = sequoia(GenoM = Geno,LifeHistData = XX, MaxSibIter = 0)
  assign('ParOUT',ParOUT, envir = globalenv())

  #Pour modifier ou non si on le desire la matrice AgePrior

  AP =  as.matrix(ParOUT$AgePriors)
  def = c("M","P","MGM","PGF","MGF","FS","MS","PS","UA")
  nb_row = nrow(AP)
  correct = FALSE
  print(AP)
  line_ap = readline('Correct AgePriors ? (y/Enter) ')
  if (line_ap != ''){
    correct = TRUE
    }
  while( line_ap != ''){
  line_ap = readline('AgePriors row to correct : '  )
  line_ap = as.numeric(line_ap)
  while(is.na(line_ap) == TRUE | (line_ap > nb_row | line_ap < 1)){
    line_ap = as.numeric(readline('Integer needed or value out of range try again. Row : '))
  }

  if (line_ap != ''){
      col_ap = readline('Column : ')
      while (col_ap %nin% def){
        col_ap = readline('Impossible value, try again. Column : ')
      }
      val_ap = as.numeric(readline('Value : '))
      while (is.na(val_ap) == TRUE){
        val_ap = as.numeric(readline('Impossible value, try again. Value : '))
      }
      AP[line_ap, col_ap] = val_ap
      assign('AP',AP, envir = globalenv())
    }
  print(AP)
  line_ap = readline('Next/Stop? (y/Enter)')
  }
  age_prior2 = readline('Used already built AgePriors ? (y/Enter) ')
  if (age_prior2 != ''){
    try({age_prior2 = readline('AgePrior file name ? ')
         AP = read.table(age_prior2, header = TRUE, sep = '\t')
         row.names(AP) = as.character(array(1:nrow(AP)))
       })
  }
  AP = as.matrix(AP)
  assign('AP',AP, envir = globalenv())
  raw_file = str_sub(raw_file,1,-5)
  ParOUT$AgePriors = AP
  assign('ParOUT2',ParOUT, envir = globalenv())
  SeqOUT <- sequoia(GenoM = Geno,SeqList = ParOUT, MaxSibIter = MS)
  assign('SeqOUT',SeqOUT, envir = globalenv())
  print('SeqOUT Printing...')
  nom = str_sub(paste(LHData),13, -1)
  write.table(SeqOUT$Pedigree, file =paste('Pedi',MS,raw_file,nom,sep = ''), sep = '\t', row.names = FALSE, quote =FALSE)
  write.table(SeqOUT$AgePriors, file =paste('AP',MS,raw_file,nom, sep = ''), sep = '\t', row.names = FALSE, quote =FALSE)
  save(SeqOUT, file =paste('SeqOUT',MS,raw_file,nom,'.RData', sep = ''))
  print('END')

  cat("\n", file = bb)
  cat(bbtime, file = bb)
  cat("\n", file = bb)
  cat(c(raw_file,genecap_file,LHData,MS),file = bb)
  closeAllConnections()
}

#'LogDB
#'
#'Allows log-in to INRA database to pick deer data needed to build the Life Hist Data.
#'
#' @param dbn Database name
#' @param hst Host server url
#' @param usr User
#' @param pw Password
#' @return Access to INRA database
#' @export
LogDB <- function(dbn,hst,usr,pw){
  log_list = c(dbn,hst,usr,pw)
  log_list2 = c('dbn','hst','usr','pw')
  for (i in 1:4){
    assign(log_list2[i], log_list[i], envir = globalenv())
  }
}

