for (i in list.files(path= "./data/", pattern = "*.ab1")){

  plot_hist(allign(get_ref(read.abif(paste0("data/",i)))))
}
