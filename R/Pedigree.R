
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

  data=read_excel(genecap_file)
  data1=subset(data,(ani_cor %nin% c("494"))) ###328 répété !
  data=droplevels(data1)
  data$ani_cor=as.character(data$ani_cor)
  data=data[,'ani_cor']
  assign('Corresp_NoDouble',data, envir = globalenv())
  drv <- dbDriver("PostgreSQL")
  ## FULL join permet de garder toutes les lignes du tableau
  con <- dbConnect(drv, dbname="db_cefs", host="pggeodb.nancy.inra.fr", user="lgervais", password="lgervais35") ## Ã  complÃ©ter

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
      dv2$mean_age[x] =as.integer(mean(c(minimum,maximum),na.rm=T))+1
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
      else{
        dv2$mean_age[x] = 'NA'
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
  # assign('LifeHisData', LifeHisData, envir = globalenv())
  write.table(LifeHistData, file =  paste('LifeHistData',nom,".txt",sep=""), sep = '\t', row.names = FALSE, quote =FALSE)
}

#' Pedigree Creator
#'
#' Creates a Pedigree taking in paramater a raw file containing a 0/1/2 genome, the link file between genetic individuals names
#' and caught individuals names.
#' It also take the Maximum Siblings number Iterations number to relatives assignment as MS1 to the Parent-Offspring assignment
#' and MS2 to the others relatives assignement.
#'
#'
#' @param raw_file The raw file name
#' @param genecap_file The genetic-catch names links file
#' @param LHData LifeHistData file
#' @param MS1 First run maximum number of iteration
#' @param MS2 Second run maximum number of iteration
#' @return The population's SNP pedigree
#' @export
#'
PediCre <- function(raw_file,genecap_file,LHData,MS1,MS2){

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

  ParOUT = sequoia(GenoM = Geno,LifeHistData = XX, MaxSibIter = MS1)
  assign('ParOUT',ParOUT, envir = globalenv())

  #Pour modifier ou non si on le desire la matrice AgePrior

  AP = ParOUT$AgePriors
  print(AP)
  line_ap = 'X'
  while( line_ap != ''){
  line_ap = readline('Modification AgePrior ligne : '  )
    if (line_ap != ''){
      line_ap = as.integer(line_ap)
      col_ap = readline('Colonne : ')
      val_ap = readline('Valeur : ')
      val_ap = as.numeric(val_ap)
      AP[line_ap, col_ap] = val_ap
    }
  print(AP)
  }
  assign('AP',AP, envir = globalenv())

  SeqOUT = sequoia(GenoM = Geno, SeqList = ParOUT, MaxSibIter = MS2)
  assign('SeqOUT',SeqOUT, envir = globalenv())
  raw_file = str_sub(raw_file,1,-5)
  write.table(SeqOUT$Pedigree, file =paste('Pedi',MS1,MS2,raw_file, sep = '-'), sep = '\t', row.names = FALSE, quote =FALSE)
}

devtools::install_github("yourusername/myfirstpackage")
