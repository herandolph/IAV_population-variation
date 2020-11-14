
library(ggplot2)
library(reshape2)
library(cobs)

##SUMMARY: Compile and load the C libraries for faster core routines.
system("R CMD SHLIB -Wall permFDR.c")
dyn.load("permFDR.so")

##SUMMARY:
permFDR=function(
##SUMMARY: df or path_to_file containing true p-values
full_data,
##SUMMARY: column of full_data where actual p_values are.
full_column_id,
##SUMMARY: df or path_to_file containing permuted p_values
perm_data,
##SUMMARY: column or columns of perm_data where permuted p_values are
perm_column_ids,
##SUMMARY: list with alternative methods to present in output
alt_methods,
##SUMMARY:Number of knots used for the splines fit (minimum 3, maximum: min(3,number_of_points_to_fit/100)). If missing the number of knots that will be used is max(5,min(20,number_of_points_to_fit/1000))
rho_o_splines_knots,
##SUMMARY:Minimum and maximum knopt sizes for the splines fit (relative to the total x-axis range (p-values) in the fit).
rho_o_splines_knots_range,
##SUMMARY:Number of iterations to be used by the splines fit.
rho_o_splines_iterations=10000,
##SUMMARY: Method for estimating p_star: if FALSE: p_star is quickly estimated as the argmin of the quotient of the local grenander p-value densities estimators; if TRUE, it is identified from a thorough sweep (computationally exhaustive).
p_star_exhaustive_sweep=FALSE,
##SUMMARY: whether the summary plot is wanted or not.
plot=TRUE,
##SUMMARY: threshold for significance assessments
significance_threshold=0.05,
##SUMMARY: report_threshold in final files and axis limit in fdrs.pdf figure.
report_threshold=1,
##SUMMARY: a folder with such name will be created to store results in working directory.
output_name="outputs_PermFDR",
##SUMMARY: whether progress messages are wanted or not
print_messages=TRUE)
{
    {
        {
            message("#############################################################################################################")
            message("##################### permFDR: an optimized method for calculating False Discovery Rates ####################")
            message("##################### that estimates the fraction of null hypothesis using permutation tests ################")
            message("#############################################################################################################")
            message("_____________________ STEP 1: Check input data ______________________________________________________________")
        }
        
        flag_debug_august=0
        
        ##1. Declare check_error_vector, perm_fdr_plot & splines_fit internal functions.
        check_scaling=function(grid,lower,greater){
            difs=pmax(lower-greater,0)
            grid_plus_one=c(grid[2:length(grid)],grid[length(grid)])
            grid_minus_one=c(grid[1],grid[1:(length(grid)-1)])
            
            bubbles=difs*(grid_plus_one-grid_minus_one)/2
            contributing_set=which(bubbles>0)
            
            for(i in 2:length(contributing_set))
            {
                bubbles[contributing_set[i]]=bubbles[contributing_set[i]-1]+bubbles[contributing_set[i]]
            }
            
            p_end=grid[which.max(bubbles)]
            max_bubble_size=bubbles[which.max(bubbles)]
            
            set=which(bubbles[1:which.max(bubbles)]==0)
            p_begin=grid[set[length(set)]+1]
            return(list(begin=p_begin,end=p_end,scaling_error=max_bubble_size))
        }
    evaluate_rho_o=function(index_p_star,grid,rho_hat,rho_hat_fitted,F_true,F_o,rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,plot,print_messages=FALSE,trend){
            
            F<-F_true
            if(index_p_star>length(F)/2)
            {
                if(trend=="none"){
                    trend="none_left"
                    #message("Evaluating p_star ",grid[index_p_star])
                }else{trend="decrease"}
                
                entries=index_p_star
                rho_hat[1:entries]=(F[index_p_star]-F[1:index_p_star])/(F_o[index_p_star]-F_o[1:index_p_star])
                
                finite_set=which(is.finite(as.numeric(as.character(rho_hat[1:entries]))))
                finite_set_grid=finite_set
                
                f_hat=splines_fit(grid[finite_set_grid],rho_hat[finite_set],rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,trend=trend,full_number=F[index_p_star],perm_number=F_o[index_p_star],print_messages=print_messages)
                rho_hat_fitted[1:length(finite_set)]<-f_hat$fitted
                
                rho_o_scaling=max(min(rho_hat_fitted[length(finite_set)],1),0)
                rm(f_hat)
            }else{
                if(trend=="none"){
                    trend="none_right"
                    #message("Evaluating p_star ",grid[index_p_star])654934466
                }else{trend="increase"}
                
                entries=(length(F)-index_p_star+1)
                rho_hat[1:entries]=(F[index_p_star]-F[index_p_star:length(F)])/(F_o[index_p_star]-F_o[index_p_star:length(F_o)])
                
                finite_set=which(is.finite(as.numeric(as.character(rho_hat[1:entries]))))
                finite_set_grid=finite_set+index_p_star-1
                
                f_hat=splines_fit(grid[finite_set_grid],rho_hat[finite_set],rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,trend=trend,full_number=F[index_p_star],perm_number=F_o[index_p_star],print_messages=print_messages)
                rho_hat_fitted[1:length(finite_set)]<-f_hat$fitted
                
                rho_o_scaling=max(min(rho_hat_fitted[1],1),0)
                rm(f_hat)
            }
            
            if(plot==TRUE)
            {
                rho_fit_table=data.frame(p_values_grid=grid[finite_set_grid],rho_hat=rho_hat[1:length(finite_set)],rho_hat_fitted=rho_hat_fitted[1:length(finite_set)])
                return(list(rho_o_scaling=rho_o_scaling,rho_fit_table=rho_fit_table))
            }else{
                return(list(rho_o_scaling=rho_o_scaling))}
        }
        estimate_rho_o=function(grid,F,f,F_o,f_o,print_messages,reference_size,p_star_exhaustive_sweep,
full_data_is_file,rho_o_splines_knots,rho_o_splines_knots_range,rho_o_splines_iterations){
            if(print_messages){	message("Estimating scaling point p*...")}
            
            ##If plot==TRUE, we'll use the Grenander estimators for producing the figure.
            if(p_star_exhaustive_sweep==FALSE | plot==TRUE)
            {
                grenander=.C("estimate_p_star",
                grid=grid,
                grid_size=as.integer(reference_size),
                F=F,
                F_o=F_o,
                f=f,
                f_o=f_o,
                F_g=F_g,
                F_o_g=F_o_g,
                f_g=f_g,
                f_o_g=f_o_g,
                index_p_star=as.integer(1),
                fitting_goodness=fitting_goodness,
                flag_grenander_choice=flag_grenander_choice)
                
                index_p_star=grenander$index_p_star
                fitting_goodness=grenander$fitting_goodness
                flag_grenander_choice=grenander$flag_grenander_choice
            }
            if(p_star_exhaustive_sweep==TRUE)
            {
                reduction_factor=(length(F))/100
                
                if(full_data_is_file){
                    p_star_indexes_array=unique(c(seq(1,reference_size,as.integer(reduction_factor)),reference_size))
                }else{
                    v_aux=vector(mode = "double", length = as.integer(reference_size/reduction_factor+1))
                    delta=as.double((max_stat-min_stat)/(reference_size)*reduction_factor)
                    grid_aux=seq(min_stat,max_stat,delta)
                    index_aux=1
                    #print(head(grid_aux))
                    for(i in 2:length(v_aux))
                    {
                        while(grid_aux[i]>grid[index_aux]) {index_aux=index_aux+1}
                        
                        v_aux[i]=index_aux
                        
                    }
                    v_aux=c(v_aux[2:length(v_aux)],reference_size)
                    p_star_indexes_array=unique(v_aux)
                }
                
                rho_o_scaling_vector=p_star_indexes_array
            
            rho_o_scaling_vector=lapply(p_star_indexes_array,evaluate_rho_o,
                grid=grid,
                rho_hat=grid,
                rho_hat_fitted=grid,
                F_true=F,
                F_o=F_o,
                rho_o_splines_knots=rho_o_splines_knots,
                rho_o_splines_iterations=rho_o_splines_iterations,
                rho_o_splines_knots_range=rho_o_splines_knots_range,
                plot=FALSE,
                print_messages=FALSE,
                trend="none")
                rho_o_scaling_vector=unlist(rho_o_scaling_vector)
                
                ##The entries of rho_o_scaling_vector are obtained from trend="none": necessary when the initial guess is far from the minimum (there, one no has a clue of how the trend could be), but anticonservative in the vicinity of it.
                ##Now, I re-calculate the lowest entries with trend=monotonic increase/decrease to get more conservative results that are in turn more comparable to what's obtained from sweep=FALSE.
                index_p_star=p_star_indexes_array[which.min(rho_o_scaling_vector)]
                flag_trend_shifted=rep(0,length(p_star_indexes_array))
                min_entry=which.min(rho_o_scaling_vector)
                cont=0
                repeat{
                    cont=cont+1
                    index_p_star=p_star_indexes_array[min_entry]
                    #print(paste(cont,min_entry,index_p_star,rho_o_scaling_vector[min_entry]))
                    
                    result=evaluate_rho_o(index_p_star,grid,rho_hat=grid,
                    rho_hat_fitted=grid,F,F_o,rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,plot=TRUE,print_messages=FALSE,trend="none")
                    
                    flag_trend_shifted[min_entry]=1
                    rho_o_scaling_vector[min_entry]=result$rho_o_scaling
                    min_entry=which.min(rho_o_scaling_vector)
                    if(flag_trend_shifted[min_entry]==1){break}
                }
                index_p_star=p_star_indexes_array[which.min(rho_o_scaling_vector)]
            }
            
            result=evaluate_rho_o(index_p_star,grid,rho_hat=grid,
            rho_hat_fitted=grid,F,F_o,rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,plot=TRUE,print_messages=print_messages,trend="undeclared")
            
            p_star=grid[index_p_star]
            rho_o_scaling=result$rho_o_scaling
            
            #Basic output
            output=list(rho_o_scaling=rho_o_scaling,index_p_star=index_p_star,p_star=p_star)
            
            if(plot==TRUE)
            {
                F_g=grenander$F_g
                F_o_g=grenander$F_o_g
                f_g=grenander$f_g
                f_o_g=grenander$f_o_g
                rho_fit_table=result$rho_fit_table
                
                output=c(output,list(flag_grenander_choice=flag_grenander_choice,fitting_goodness=fitting_goodness,F_g=F_g,F_o_g=F_o_g,f_g=f_g,f_o_g=f_o_g,rho_fit_table=rho_fit_table))
                
                if(p_star_exhaustive_sweep==TRUE)
                {
                    rho_o_scaling_loop=data.frame(p_star_candidate=grid[p_star_indexes_array],rho_o_scaling_candidate=rho_o_scaling_vector)
                    output=c(output,list(rho_o_scaling_loop_table=rho_o_scaling_loop))
                }
            }else if(p_star_exhaustive_sweep==FALSE)
            output=c(output,list(flag_grenander_choice=flag_grenander_choice,fitting_goodness=fitting_goodness))
            return(output)
        }
        perm_fdr_plot=function(plots){
            
            rel_heights <- c(1,1,1,1)
            rel_widths <- 1
            scale <- c(1,1,1,1)
            rows = 4
            cols = 1
            
            grobs <- lapply(plots, function(x) {ggplot2::ggplotGrob(x)})
            num_widths <- unique(lapply(grobs, function(x) {
                length(x$widths)
            }))
            num_widths[num_widths == 0] <- NULL
            max_widths <- do.call(grid::unit.pmax, lapply(grobs,
            function(x) {
                x$widths
            }))
            
            for (i in 1:4) {
                grobs[[i]]$widths <- max_widths
            }
           
            x_deltas <- rel_widths/sum(rel_widths)
            y_deltas <- rel_heights/sum(rel_heights)
            xs <- cumsum(rel_widths)/sum(rel_widths) - x_deltas
            ys <- 1 - cumsum(rel_heights)/sum(rel_heights)
            ggdraw=function (plot = NULL)
            {
                theme_nothing=function(base_size = 12, base_family = "")
                {
                    theme_grey(base_size = 12, base_family = "") %+replace%
                    theme(rect = element_rect(fill = "transparent", colour = NA,
                    color = NA, size = 0, linetype = 0), line = element_blank(),
                    text = element_blank(), title = element_blank(),
                    panel.background = element_blank(), axis.ticks.length = grid::unit(0,
                    "lines"), legend.position = "none", plot.margin = grid::unit(c(0,
                    0, 0, 0), "lines"))
                }
                
                d <- data.frame(x = 0:1, y = 0:1)
                p <- ggplot(d, aes_string(x = "x", y = "y")) + scale_x_continuous(limits = c(0,
                1), expand = c(0, 0)) + scale_y_continuous(limits = c(0,
                1), expand = c(0, 0)) + theme_nothing() + labs(x = NULL,y = NULL)
                if (!is.null(plot)) {
                    if (methods::is(plot, "ggplot")){
                        g<-ggplot2::ggplotGrob(plot)
                    }else if (methods::is(plot, "gtable")){
                        g<-plot
                    }else{
                        stop('Argument needs to be of class "ggplot" or "gtable"' )
                    }
                    plot.grob <- grid::grobTree(g)
                    p <- p + annotation_custom(plot.grob)
                }
                p
            }
            p <- ggdraw()
            col_count <- 0
            row_count <- 1
            for (i in 1:(rows * cols)) {
                if (i > 4)
                break
                x_delta <- x_deltas[col_count + 1]
                y_delta <- y_deltas[row_count]
                width <- x_delta * scale[i]
                height <- y_delta * scale[i]
                x_off <- (x_delta - width)/2
                y_off <- (y_delta - height)/2
                x <- xs[col_count + 1] + x_off
                y <- ys[row_count] + y_off
                draw_grob=function (grob, x = 0, y = 0, width = 1, height = 1)
                {
                    annotation_custom(grid::grobTree(grob), xmin = x, xmax = x + width, ymin = y, ymax = y + height)
                }
                
                p <- p + draw_grob(grid::grobTree(grobs[[i]]),x,y,width,height)
                col_count <- col_count + 1
                if (col_count >= cols) {
                    col_count <- 0
                    row_count <- row_count + 1
                }
            }
            p
        }
    splines_fit=function(stats_array,rho_hat,rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,trend,full_number,perm_number,print_messages=FALSE,constraint_matrix=NULL){
           
           if(rho_o_splines_knots==0){
               #rho_o_splines_knots=min(30,as.integer(2+(28/3)*log(length(stats_array)/10)))
               rho_o_splines_knots=min(20,as.integer(2+3*log10(length(stats_array)/10)))
               if(print_messages==TRUE)message("Fitting rho_o from a spline with ",rho_o_splines_knots," knots...")
           }
           if(all(rho_o_splines_knots_range==c(0,0))){
               #rho_o_splines_knots_range=c(1e-12,1)
               rho_o_splines_knots_range=c(0,1)
               
           }
           
           p_star=stats_array[length(stats_array)]
           p_star_index=length(stats_array)
           
           if(trend=="decrease"|trend=="none_left")
           {
               if(missing(constraint_matrix))
               constraint_matrix=as.matrix(data.frame(c(2),c(p_star),c(0)))
               width_last_one_per_thousand_quantile=p_star-stats_array[as.integer(length(stats_array)*0.999)]
           }else if(trend=="increase"|trend=="none_right"){
               if(missing(constraint_matrix))
               constraint_matrix=as.matrix(data.frame(c(2),c(stats_array[1]),c(0)))
               width_last_one_per_thousand_quantile=stats_array[as.integer(length(stats_array)*0.001)]-stats_array[1]
           }
           
           if(trend %in% c("none_left","none_right")) trend="none"
           if(trend %in% c("decrease","increase")) trend="convex"


           df<-data.frame(stats_array,rho_hat)
           #write.table(df,paste0(output_name,"/cobs_input.txt"))
           
           range=abs(stats_array[1]-stats_array[length(stats_array)])
           max_res=range*max(rho_o_splines_knots_range)
           min_res=range*min(rho_o_splines_knots_range)
           
           if((trend =="convex") & max_res<width_last_one_per_thousand_quantile)
           {
               max_res=width_last_one_per_thousand_quantile
               if(print_messages==TRUE)message(paste0("(Note: maximum knot size augmented to ",max_res," from proposed input to allow the closest knot to p_star to cover the 1 per 1000 of the data at least)"))
           }
           nknots=rho_o_splines_knots
           
           repeat{
               knots=stats_array[1]
               index_last_knot=1
               step=as.integer(p_star_index/nknots)
               
               repeat
               {
                   past_index=index_last_knot
                   keep_going=(index_last_knot+step)<p_star_index
                   if(keep_going==1){
                       anchura_quantile=stats_array[as.integer(index_last_knot+step)]-knots[length(knots)]
                       end=stats_array[index_last_knot+step]
                       end_index=index_last_knot+step
                   }else{
                       anchura_quantile=p_star-knots[length(knots)]
                       end=p_star
                       end_index=p_star_index
                   }
                   if(anchura_quantile<min_res){
                       knots=c(knots,min(knots[length(knots)]+min_res,p_star))
                       if(length(which(stats_array>=knots[length(knots)]))>0)
                       {
                           index_last_knot=which(stats_array>=knots[length(knots)])[1]
                       }else{
                           index_last_knot=p_star_index
                           break;
                       }
                   }else{
                       if(anchura_quantile>max_res){
                           knots=c(knots,knots[length(knots)]+max_res)
                           if(length(which(stats_array>=knots[length(knots)]))>0)
                           {
                               index_last_knot=which(stats_array>=knots[length(knots)])[1]
                           }else{
                               index_last_knot=p_star_index
                               break;
                           }
                       }else{
                           knots<-c(knots,end)
                           index_last_knot=end_index
                           if(end==1){break}
                       }
                   }
                   if((p_star-knots[length(knots)])<0.2*(knots[length(knots)]-knots[length(knots)-1]))
                   {
                       knots[length(knots)]=p_star
                       index_last_knot=length(knots)
                       break
                   }
               }
               if(length(knots)>=rho_o_splines_knots){
                   break
               }else{nknots=nknots+1}
           }
           
            f_hat<-suppressWarnings(cobs(stats_array,rho_hat,knots=knots,constraint=trend,pointwise=constraint_matrix,maxiter=rho_o_splines_iterations,print.warn=FALSE,print.mesg=FALSE))
           return(f_hat)
       }
       
        check_error_vector=function(error_vector,flag_error){
            for(i in 1:length(error_vector))
            {
                if(nchar(error_vector[i])>0)
                {
                    message(error_vector[i])
                    flag_error=1
                }
            }
        }
        
    
        ##SUMMARY 2. Check plot output_name format: possible warning: it's not a single logical value: will be set to TRUE then.
        {
            warning_output_name=vector(mode = "character", length = 1)
            if(typeof(output_name)!="character"| length(output_name)!=1)
            {
                warning_output_name[1]=paste0("WARNING: Invalid argument output_name: Outputs will be saved in: ",getwd(),"/outputs_PermFDR")
                output_name="outputs_PermFDR"
            }
            if(output_name==""){
                command="mkdir -p outputs_PermFDR"
                output_name="outputs_PermFDR"
                warning_output_name[1]=paste0("WARNING: Invalid argument output_name: Outputs will be saved in: ",getwd(),"/outputs_PermFDR")
                system(command)
            }else{
                command=paste0("mkdir -p ",output_name)
                system(command)
            }
            if(nchar(warning_output_name)>0) message(warning_output_name)
        }
        ##SUMMARY: 3. Check full_data & full_column_id inputs formats:
        ##For full_data: possible errors:
        ##1. Missing argument.
        ##2. Bad-format argument (neither df nor file).
        ##3. If file, and the file is absent or impossible to load: invalid file.
        ##4: if file: bad parsing.
        ##5: Not numeric, or values outside interval [0,1] in tests column (invalid p-values).
        ##For full_column_id: possible errors.
        ##1. Missing argument.
        ##2. Bad-format argument (neither column name or ordinal).
        ##3. Out-of-range argument(column name or ordinal that doesn't match with present columns ids in full_data
        {
            error_full_data=vector(mode = "character", length = 1)
            error_full_column_id=vector(mode = "character", length = 1)
            if(missing(full_data)){
                error_full_data[1]="ERROR: Absent argument full_data."
            }else if((typeof(full_data)!="character"|(typeof(full_data)=="character" & length(full_data)>1)) & is.data.frame(full_data)==FALSE){
                error_full_data[1]="ERROR: Invalid argument full_data: must be a data frame or a path to a file."
            }else if(typeof(full_data)=="character"){
                full_data_is_file=1
                #Check wether a rownames column is present: NOTE that a colnames row is needed
                head=try(read.table(full_data,nrows=2,header=TRUE),silent=TRUE)
                head_2=try(read.table(full_data,nrows=2,header=FALSE),silent=TRUE)
                if(is(head,"try-error"))
                {error_full_data[1]=paste0("ERROR: Invalid argument full_data: file: ",full_data," is missing or bad formatted.")
                }else{
                    if(is(head_2,"try-error")) {
                        rownames_flag_full=1
                    }else{
                        rownames_flag_full=0
                    }
                    columns_full_data=ncol(head)+rownames_flag_full
                }
                
                #Check if full_column_id argument matches a column in full_data file
                if(nchar(error_full_data[1])==0)
                {
                    if(missing(full_column_id)){
                        error_full_column_id[1]="ERROR: Absent argument full_column_id."
                    }else if((typeof(full_column_id)!="character" & typeof(full_column_id)!="double") | length(full_column_id)>1){
                        error_full_column_id[1]="ERROR: Invalid argument full_column_id: must be a valid column identificator (column name or integer ordinal index)"
                    }else if(typeof(full_column_id)=="character"){
                        buffer=which(colnames(head) %in% full_column_id)
                        if(length(buffer)!=1) error_full_column_id[1]=paste0("ERROR: Invalid argument full_column_id:\n","Column name introduced: ",full_column_id," was not found in full_data file: ",full_data)
                    }else if(full_column_id%%1!=0){
                        error_full_column_id[1]="ERROR: Invalid argument full_column_id: must be a valid column identificator (column name or integer ordinal index)"
                    }else if(full_column_id<1 | full_column_id>ncol(head))
                    {
                        error_full_column_id[1]=paste0("ERROR:  Invalid argument full_column_id:\nColumn ordinal: ",full_column_id," is out of range: full_data file contains ",ncol(head)," columns.")
                    }else{
                        buffer=full_column_id
                    }
                    full_column_id=buffer+rownames_flag_full
                }
                
                #If everything is yet, now check file parsing and get number of rows:
                if(nchar(error_full_data[1])==0 & nchar(error_full_column_id[1])==0)
                {
                    fields=ncol(head)
                    if(rownames_flag_full==1) fields=fields+1
                    if(print_messages){message("Checking full_data file...")}
                    ##This will do the following:
                    ##Keep track of max_p and min_p
                    ##Print one line per bad-parsed line.
                    ## IF any value outside interval [0,1] is found standing for a p-value: flag_badp=1
                    ## Print a final 4 lines with NR (number of rows in file (number of tests plus one for the colnames)) min_p, max_p, flag_badp
                    
                    command=paste0("awk 'BEGIN{RS=\"\\r\";print 1;}{if(NR>1){print NR; exit 0;}}' ",full_data,"> ",output_name,"parsing_full1.txt")
                    system(command)
                    full_parsing_name<-paste0(output_name,"parsing_full1.txt")
                    parsing=read.table(full_parsing_name,header=FALSE)
                    system(paste0("rm ",full_parsing_name))
                    if(nrow(parsing)==2){
                        command=paste0("awk 'BEGIN {max=",fields,";flag_bad_p=0;RS=\"\\r\";min_stat=1;max_stat=0} {if(NF!=max) {print NR;}if(NR>1){if(($",full_column_id,"<0)||($",full_column_id,"!~ '/^[0-9.eE+'-']*$/')||($",full_column_id,">1))  {flag_bad_p=1; exit 0;}} if(NR>1) {if(($",full_column_id,"<min_stat)&&($",full_column_id,">0))  min_stat=$",full_column_id,";if($",full_column_id,">max_stat) max_stat=$",full_column_id,";}} END{print NR; print min_stat; print max_stat; print flag_bad_p}' ",full_data,"> ",output_name,"parsing_full.txt\n")
                        system(command)
                    }else{
                        command=paste0("awk 'BEGIN{RS=\"\\n\";print 1;}{if(NR>1){print NR; exit 0;}}' ",full_data,"> ",output_name,"parsing_full2.txt")
                        system(command)
                        full_parsing_name<-paste0(output_name,"parsing_full2.txt")
                        parsing=read.table(full_parsing_name,header=FALSE)
                        system(paste0("rm ",full_parsing_name))
                        if(nrow(parsing)==2){
                           command=paste0("awk 'BEGIN {max=",fields,";flag_bad_p=0;RS=\"\\n\";min_stat=1;max_stat=0} {if(NF!=max) {print NR;}if(NR>1){if(($",full_column_id,"<0)||($",full_column_id,"!~ '/^[0-9.eE+'-']*$/')||($",full_column_id,">1))  {flag_bad_p=1; exit 0;}} if(NR>1) {if(($",full_column_id,"<min_stat)&&($",full_column_id,">0))  min_stat=$",full_column_id,";if($",full_column_id,">max_stat) max_stat=$",full_column_id,";}} END{print NR; print min_stat; print max_stat; print flag_bad_p}' ",full_data,"> ",output_name,"parsing_full.txt\n")
                           system(command)
                        }else{
                            error_full_data[1]="ERROR: Invalid argument full_data"
                        }
                    }
                    full_parsing_name<-paste0(output_name,"parsing_full.txt")
                    parsing=read.table(full_parsing_name,header=FALSE)
                    system(paste0("rm ",full_parsing_name))
                    ##Lines that the file will have by default: the final 4 lines plus another one (indicating that the first row contains one field less than the default) if a rownames column is present.
                    parsing_lines=4+rownames_flag_full
                    if(nrow(parsing)>parsing_lines)
                    {error_full_data[1]="ERROR: Invalid argument full_data: bad parsing."
                    }else{
                        full_tests=as.double(parsing[nrow(parsing)-3,1])-1
                        
                        min_stat=parsing[nrow(parsing)-2,1]
                        max_stat=parsing[nrow(parsing)-1,1]
                        if(parsing[nrow(parsing),1]==1)
                        {error_full_data[1]=paste0("ERROR: Invalid p-values (non-numeric and/or out of range) at full_data file:",full_data)}
                    }
                    rm(parsing)
                }
            }else if (is.data.frame(full_data)==TRUE){
                
                full_data_is_file=0
                
                #Check if full_column_id argument matches a column in full_data dataframe
                if(missing(full_column_id)){
                    error_full_column_id[1]="ERROR: Absent argument full_column_id."
                }else if((typeof(full_column_id)!="character" & typeof(full_column_id)!="double") | length(full_column_id)>1){
                    error_full_column_id[1]="ERROR: Invalid argument full_column_id: must be a valid column identificator (column name or integer ordinal index)"
                }else if(typeof(full_column_id)=="character")
                {
                    buffer=which(colnames(full_data) %in% full_column_id)
                    if(length(buffer)!=1)
                    error_full_column_id[1]=paste0("ERROR: Invalid argument full_column_id:\nColumn name introduced: ",full_column_id," was not found in full_data dataframe.")
                }else if(full_column_id%%1!=0){
                    error_full_column_id[1]="ERROR: Invalid argument full_column_id: must be a valid column identificator (column name or integer ordinal index)"
                }else if(full_column_id<1 | full_column_id>ncol(full_data))
                {
                    error_full_column_id[1]=paste0("ERROR:  Invalid argument full_column_id:\nColumn ordinal: ",full_column_id," is out of range: full_data file contains ",ncol(full_data)," columns.")
                }else{
                    buffer=full_column_id
                }
                
                if(nchar(error_full_column_id[1])==0)
                {
                    full_column_id=buffer
                    failed_full_entries=length(which(is.na(full_data[,full_column_id]) | (full_data[,full_column_id]<0) | (full_data[,full_column_id]>1)))
                    if(failed_full_entries>0) {error_full_data[1]="ERROR: Invalid p-values at full_data dataframe (non-numeric and/or out of range)."}
                }
                
                if(nchar(error_full_column_id[1])==0 & nchar(error_full_data[1])==0)
                {
                    full_tests=nrow(full_data)
                }
            }
            check_error_vector(error_full_data,flag_error)
            check_error_vector(error_full_column_id,flag_error)
        }
        ##SUMMARY: 4 Check perm_data & perm_column_ids inputs formats:
        ##For perm_data: possible errors:
        ##1. Missing argument.
        ##2. Bad-format argument (neither df nor file).
        ##3. If file, and the file is absent or impossible to load: invalid file.
        ##4: if file: bad parsing.
        ##5: Not numeric, or values outside interval [0,1] in tests column (invalid p-values).
        ##For perm_column_ids: possible errors.
        ##1. Missing argument.
        ##2. Bad-format argument (neither vector of column names or ordinals or atring "all").
        ##3. The argument contains out-of-range entries (column name or ordinal that doesn't match with present columns ids in perm_data
        {
            error_perm_data=vector(mode = "character", length = 1)
            error_perm_column_id=vector(mode = "character", length = 1)
            if(missing(perm_data)){
                error_perm_data[1]="ERROR: Absent argument perm_data."
            }else if((typeof(perm_data)!="character"|(typeof(perm_data)=="character" & length(perm_data)>1)) & is.data.frame(perm_data)==FALSE){
                error_perm_data[1]="ERROR: Invalid argument perm_data: must be a data frame or a path to a file."
            }else if(typeof(perm_data)=="character"){
                perm_data_is_file=1
                #Check wether a rownames column is present: NOTE that a colnames row is needed
                head=try(read.table(perm_data,nrows=2,header=TRUE),silent=TRUE)
                head_2=try(read.table(perm_data,nrows=2,header=FALSE),silent=TRUE)
                if(is(head,"try-error"))
                {error_perm_data[1]=paste0("ERROR: Invalid argument perm_data: file: ",perm_data,"missing file or bad formatted.")
                }else{
                    if(is(head_2,"try-error")){
                        rownames_flag_perm=1
                    }else{
                        rownames_flag_perm=0
                    }
                    columns_perm_data=ncol(head)+rownames_flag_perm
                }
                
                #Check if all entries in perm_column_id argument match columns in perm_data file
                if(nchar(error_perm_data[1])==0)
                {
                    if(missing(perm_column_ids)){
                        error_perm_column_ids[1]="ERROR: Absent argument perm_column_ids."
                    }else if(typeof(perm_column_ids)!="character" & (typeof(perm_column_ids)!="double")){
                        error_perm_column_ids[1]="ERROR: Invalid argument perm_column_ids: must be a vector of valid column identificators (column names, integer ordinal indexes or \"all\")"
                    }else if(typeof(perm_column_ids)=="character")
                    {
                        if(perm_column_ids=="all")
                        {buffer=c(1:(columns_perm_data-rownames_flag_perm))
                        }else{
                            buffer=which(colnames(head) %in% perm_column_ids)
                            not_found=which(!(perm_column_ids %in% colnames(head)))
                            if(length(not_found)>0)
                            error_perm_column_id[1]="ERROR: Invalid entries at perm_column_ids: names of columns introduced not found in perm_data."
                        }
                    }else if(sum(perm_column_ids%%1)!=0){
                        error_perm_column_id[1]="ERROR: Invalid argument perm_column_ids: must be a vector of valid column identificators (column names, integer ordinal indexes or \"all\")"
                    }else{
                        valid_set=which((perm_column_ids>0)&(perm_column_ids<(ncol(head)+1)))
                        buffer=perm_column_ids[valid_set]
                        if(length(valid_set)<length(perm_column_ids))
                        error_perm_column_id[1]="ERROR: Invalid entries at perm_column_ids: out of range column ordinals introduced."
                    }
                    perm_column_ids=buffer+rownames_flag_perm
                }
                
                #Check file parsing and get number of rows:
                if(nchar(error_perm_data[1])==0 & nchar(error_perm_column_id[1])==0)
                {
                    fields=ncol(head)
                    if(rownames_flag_perm==1) fields=fields+1
                    if(print_messages){message("Checking perm_data file...")}
                    
                    command=paste0("awk 'BEGIN{RS=\"\\r\";print 1;}{if(NR>1){print NR; exit 0;}}' ",perm_data,"> ",output_name,"parsing_perm1.txt\n")
                    system(command)
                    perm_parsing_name<-paste0(output_name,"parsing_perm1.txt")
                    parsing=read.table(perm_parsing_name,header=FALSE)
                    system(paste0("rm ",perm_parsing_name))

                    if(nrow(parsing)==2){
                        command=paste0("awk 'BEGIN {max=",fields,";flag_bad_p=0;RS=\"\\r\";} {if(NF!=max) {print NR;}if(NR>1){if(($",perm_column_ids[1],"<0)||($",perm_column_ids[1],"!~ '/^[0-9.eE+'-']*$/')||($",perm_column_ids[1],">1)")
                        if(length(perm_column_ids)>1){
                            for(i in 2:length(perm_column_ids))
                            command=paste0(command,"||($",perm_column_ids[i],"<0)||($",perm_column_ids[i],"!~ '/^[0-9.e+'-']*$/')||($",perm_column_ids[i],">1) ")
                        }
                        command=paste0(command,"){flag_bad_p=1; exit 0;}}} END{print NR; print flag_bad_p}' ",perm_data,"> ",output_name,"parsing_perm.txt\n")
                        system(command)
                    }else{
                        command=paste0("awk 'BEGIN{RS=\"\\n\";print 1;}{if(NR>1){print NR; exit 0;}}' ",perm_data,"> ",output_name,"parsing_perm2.txt\n")
                        system(command)
                        perm_parsing_name<-paste0(output_name,"parsing_perm2.txt")
                        parsing=read.table(perm_parsing_name,header=FALSE)
                        system(paste0("rm ",perm_parsing_name))
                        if(nrow(parsing)==2){
                            command=paste0("awk 'BEGIN {max=",fields,";flag_bad_p=0;RS=\"\\n\";} {if(NF!=max) {print NR;}if(NR>1){if(($",perm_column_ids[1],"<0)||($",perm_column_ids[1],"!~ '/^[0-9.eE+'-']*$/')||($",perm_column_ids[1],">1)")
                            if(length(perm_column_ids)>1){
                                for(i in 2:length(perm_column_ids))
                                command=paste0(command,"||($",perm_column_ids[i],"<0)||($",perm_column_ids[i],"!~ '/^[0-9.e+'-']*$/')||($",perm_column_ids[i],">1) ")
                            }
                            command=paste0(command,"){flag_bad_p=1; exit 0;}}} END{print NR; print flag_bad_p}' ",perm_data,"> ",output_name,"parsing_perm.txt\n")
                            system(command)
                        }else{
                            error_perm_data[1]="ERROR: Invalid argument perm_data"
                        }
                    }
                    
                    perm_parsing_name<-paste0(output_name,"parsing_perm.txt")
                    parsing=read.table(perm_parsing_name,header=FALSE)
                    system(paste0("rm ",perm_parsing_name))
                    parsing_lines=2+rownames_flag_perm
                    if(nrow(parsing)>parsing_lines)
                    {
                        error_perm_data[1]="ERROR: Invalid argument perm_data: bad parsing."
                        return(0)
                    }else
                    {
                        perm_tests=(parsing[nrow(parsing)-1,1]-1)*length(buffer)
                        
                        if(parsing[nrow(parsing),1]==1)
                        {error_perm_data[1]="ERROR: Invalid p-values at perm_data file (non-numeric and/or out of range)."}
                        
                        rm(parsing)
                    }
                }
            } else if (is.data.frame(perm_data)==TRUE){
                perm_data_is_file=0
                #Check if all entries in perm_column_id argument match columns in perm_data dataframe
                if(missing(perm_column_ids)){
                    error_perm_column_ids[1]="ERROR: Absent argument perm_column_ids."
                }else if(typeof(perm_column_ids)!="character" & (typeof(perm_column_ids)!="double")){
                    error_perm_column_ids[1]="ERROR: Invalid argument perm_column_ids: must be a vector of valid column identificators (column names, integer ordinal indexes or \"all\")"
                }else if(typeof(perm_column_ids)=="character")
                {
                    if(perm_column_ids=="all")
                    {
                        buffer=c(1:ncol(perm_data))
                    }else{
                        buffer=which(colnames(perm_data) %in% perm_column_ids)
                        not_found=which(!(perm_column_ids %in% colnames(perm_data)))
                        if(length(buffer)==0 | length(not_found)>0)
                        error_perm_column_id[1]="ERROR: Invalid entries at perm_column_ids: names of columns introduced not found in perm_data."
                    }
                    
                }else if(sum(perm_column_ids%%1)!=0){
                    error_perm_column_id[1]="ERROR: Invalid argument perm_column_ids: must be a vector of valid column identificators (column names, integer ordinal indexes or \"all\")"
                }else
                {
                    valid_set=which((perm_column_ids>0)&(perm_column_ids<(ncol(perm_data)+1)))
                    buffer=perm_column_ids[valid_set]
                    if(length(valid_set)<length(perm_column_ids))
                    error_perm_column_id[1]="ERROR: Invalid entries at perm_column_ids:out of range column ordinals introduced."
                }
                perm_column_ids=buffer
                
                failed_perm_entries=length(which(is.na(perm_data[,perm_column_ids]) | (perm_data[,perm_column_ids]<0) | (perm_data[,perm_column_ids]>1)))
                if(failed_perm_entries>0) {error_perm_data[1]="ERROR: Invalid p-values at full data file (non-numeric and/or out of range)."}
                
                if(nchar(error_perm_column_id[1])==0 & nchar(error_perm_data[1])==0)
                {
                    #when a df is provided all data is sampled:
                    #load only columns of interest in perm_data into a vector
                    perm_data=unlist(perm_data[,perm_column_ids])
                    perm_tests=length(perm_data)
                    perm_data=perm_data[order(as.double(perm_data))]
                    perm_data<-as.double(as.character(perm_data))
                }
            }
            check_error_vector(error_perm_data,flag_error)
            check_error_vector(error_perm_column_id,flag_error)
        }
        
        ##SUMMARY: 5 Check alt_methods input format: possible warning: bad format, none or some of the entries do not match one of the possible values: alt_methods set to the subset of valid entries, if any.
        {
            warning_alt_methods=vector(mode = "character", length = 1)
            valid_alt_methods=c("Fdr_ST","Fndr_ST","Fdr_BH_unif","Fdr_BH_perm","Fdr_MV","BF")
            if(missing(alt_methods)){alt_methods="none"}
            set_alt_methods=which(valid_alt_methods %in% alt_methods)
            if(typeof(alt_methods)!="character")
            {
                alt_methods="none"
                warning_alt_methods[1]="WARNING: Invalid argument alt_methods: if any, it must be a vector containing some of the alternative methods ids: 'Fdr_ST','Fndr_ST','Fdr_BH_unif','Fdr_BH_perm','Fdr_MV','BF'; no alternative fdrs estimator will be computed."
            }else if(length(alt_methods)>length(set_alt_methods) & alt_methods!="none"){
                if(length(set_alt_methods)==0){
                    alt_methods="none"
                    warning_alt_methods[1]="WARNING: Invalid entries at alt_methods: alt_methods set to 'none'"
                }else{
                    alt_methods=valid_alt_methods[set_alt_methods]
                    warning_alt_methods[1]="WARNING: Invalid entries at alt_methods; alt_methods set to:\n"
                    for(i in 1:length(alt_methods))
                    warning_alt_methods[1]=paste0(warning_alt_methods[1],"'",alt_methods[i],"' ")
                    warning_alt_methods[1]=paste0(warning_alt_methods[1],".")
                }
            }
            if(nchar(warning_alt_methods)>0) message(warning_alt_methods)
        }
        ##SUMMARY: 6 Check pi_o_knot and rho_o_splines_knots_range formats: possible warnings:
        ## rho_o_splines_knots it's not a single integer value: will be set to 0 then (to be chose authomatically)
        ##And or rho_o_splines_knots_range is not a vector of two numbers between 0 and 1: set to zero then (to be chose authomatically)
        
        warning_rho_o_splines_knots=vector(mode = "character", length = 1)
        
        if(missing(rho_o_splines_knots))
        {
            rho_o_splines_knots=0
        }else if(!is.numeric(rho_o_splines_knots)|length(rho_o_splines_knots)!=1|rho_o_splines_knots%%1!=0|rho_o_splines_knots<2|rho_o_splines_knots>full_tests/10){
            warning_rho_o_splines_knots[1]=paste0("WARNING: invalid argument rho_o_splines_knots: must be a single integer value between 2 and ",full_tests/10,".\nThe number of knots will be calculated automatically")
            rho_o_splines_knots=0
        }
        
        if(!missing(rho_o_splines_knots_range)){
            if(!is.numeric(rho_o_splines_knots_range)| length(rho_o_splines_knots_range)!=2 | length(which(rho_o_splines_knots_range<0|rho_o_splines_knots_range>1))>0){
                if(length(warning_rho_o_splines_knots)>1)
                {
                    warning_rho_o_splines_knots[1]="WARNING: Invalid pi_o_knot and rho_o_splines_knots_range arguments\nthe number and sizes of knots for rho_o fit will be chosen automatically."
                }else{
                    warning_rho_o_splines_knots[1]="WARNING: Invalid rho_o_splines_knots_range argument\nthe allowed knot range for rho_o fit will be chosen automatically."}
                rho_o_splines_knots_range=c(0,0)
            }
        }else{rho_o_splines_knots_range=c(0,0)}
        
        if(nchar(warning_rho_o_splines_knots)>0) message(warning_rho_o_splines_knots)
        
        ##SUMMARY: 7 Check rho_o_splines_iterations input format: possible warning: it's not a single numeric value, or it is not integer. In that case, put to 10000.
        {
            warning_rho_o_splines_iterations=vector(mode = "character", length = 1)
            if(!is.numeric(rho_o_splines_iterations)|length(rho_o_splines_iterations)!=1|rho_o_splines_iterations%%1!=0)
            {
                warning_rho_o_splines_iterations[1]=paste0("WARNING: invalid argument rho_o_splines_iterations: must be a single integer value. Set to 10000")
                rho_o_splines_iterations=10000
            }
            if(nchar(warning_rho_o_splines_iterations)>0) message(warning_rho_o_splines_iterations)
        }
        ##SUMMARY: 8 Check p_star_exhaustive_sweep input format: possible warning: it's not a single logical value: will be set to TRUE then.
        {
            warning_p_star_exhaustive_sweep=vector(mode = "character", length = 1)
            if(typeof(p_star_exhaustive_sweep)!="logical"| length(p_star_exhaustive_sweep)!=1)
            {
                warning_p_star_exhaustive_sweep[1]=paste0("WARNING: invalid argument p_star_exhaustive_sweep: must be a logical variable, set to FALSE.")
                p_star_exhaustive_sweep=FALSE
            }
            if(nchar(warning_p_star_exhaustive_sweep)>0) message(warning_p_star_exhaustive_sweep)
        }
        ##SUMMARY: 9 Check plot input format: possible warning: it's not a single logical value: will be set to TRUE then.
        {
            warning_plot=vector(mode = "character", length = 1)
            if(typeof(plot)!="logical"| length(plot)!=1)
            {
                warning_plot[1]=paste0("WARNING: invalid argument plot: must be a logical variable, set to TRUE.")
                plot=TRUE
            }
            if(nchar(warning_plot)>0) message(warning_plot)
        }
        ##SUMMARY: 10 Check significance_threshold input format: possible warning: it's not a single numeric value value or it is not within (0,1]. Set to 0.05 then
        {
            warning_significance_threshold=vector(mode = "character", length = 1)
            if(!is.numeric(significance_threshold)|length(significance_threshold)!=1)
            {
                warning_significance_threshold[1]=paste0("ERROR: invalid argument significance_threshold: must be a single numeric value in (0,1]. Set to 0.05")
                significance_threshold=0.05
            }else if(significance_threshold<=0 |significance_threshold>1){
                warning_significance_threshold[1]=paste0("ERROR: invalid argument significance_threshold: out of range: it must be a numeric value in (0,1]. Set to 0.05")
                significance_threshold=0.05
            }
            if(nchar(warning_significance_threshold)>0) message(warning_significance_threshold)
        }
        ##SUMMARY: 11 Check report_threshold input format: possible warning: it's not a single numeric value value or it is not within (0,1]. Set to 1 then.
        {
            warning_report_threshold=vector(mode = "character", length = 1)
            if(!is.numeric(report_threshold)|length(report_threshold)!=1)
            {
                warning_report_threshold[1]=paste0("ERROR: invalid argument report_threshold: must be a single numeric value in (0,1]. Set to 1")
                report_threshold=1
            }else if(report_threshold<=0 |report_threshold>1){
                warning_report_threshold[1]=paste0("ERROR: invalid argument report_threshold: out of range: it must be a numeric value in (0,1]. Set to 1")
                report_threshold=1
            }
            if(nchar(warning_report_threshold)>0) message(warning_report_threshold)
        }
        ##SUMMARY: 12 Check print_messages input format: possible warning: it's not a single logical value: will be set to TRUE then.
        {
            warning_print_messages=vector(mode = "character", length = 1)
            if(typeof(print_messages)!="logical"| length(print_messages)!=1)
            {
                warning_output_name[1]=paste0("WARNING: invalid argument print_messages: must be a logical variable, set to TRUE.")
                print_messages=TRUE
            }
            if(nchar(warning_print_messages)>0) message(warning_print_messages)
        }
        
        if(exists("flag_error")) return(0)
        
        #Other warnings:
        warning_full_size=vector(mode = "character", length = 1)
        warning_perm_size=vector(mode = "character", length = 1)
        
        ##There are less than 1000 tests (pi_o estimation loses roboustness, makes little sense to run the algorithm on files)
        if(typeof(full_data)=="character")
        {
            if(full_tests<1000)
            {
                if(plot==FALSE)
                {
                    plot=TRUE
                    warning_full_size[1]=paste0("WARNING: few tests, you might want to load full_data from data.frame for non-interpolated fdr calculations.\nWARNING: pi_o stimations are more robust for larger number of tests. Check distributions plots. (plot set to TRUE).")
                }else{
                    warning_full_size[1]=paste0("WARNING: few tests, you might want to load full_data from data.frame for non-interpolated fdr calculations.\nWARNING: pi_o stimations are more robust for larger number of tests. Check distributions plots.")
                }
            }else if(full_tests<100000){
                {
                    warning_full_size[1]=paste0("WARNING: moderate number of tests: you might want to load full_data from data.frame for non-interpolated fdr calculations.")
                }
            }
        }else if(full_tests<1000)
        {
            if(plot==FALSE)
            {
                plot=TRUE
                warning_full_size[1]=paste0("WARNING: pi_o stimations are more robust for larger number of tests. Check distributions plots (plot set to TRUE).")
            }else{
                warning_full_size[1]=paste0("WARNING: few tests, pi_o stimations are more robust for larger number of tests. Check distributions plots.")
            }
        }
        ##There are less than 10 permutations of the true data (might not be enough: fdr_Bootstrap recommended)
        if(perm_tests<10*full_tests)
        {
            warning_perm_size[1]=paste0("WARNING: reduced number of permutations: increasing the number of permutation tests might result in more robust results (Try fdr_Bootstrap)")
        }
        if(nchar(warning_full_size)>0) message(warning_full_size)
        if(nchar(warning_perm_size)>0) message(warning_perm_size)
    }
    if(print_messages)
    {
        message("\n_____________________ STEP 2: Obtain p-value distributions distibutions _____________________________________")
    }
    if(print_messages)
    {    message("Obtaining experimental p-values distribution...") }
    if(full_data_is_file)
    {
        reference_size=min(full_tests,100000)

        grid_unif=vector(mode = "double", length = reference_size)
        grid_log=vector(mode = "double", length = reference_size)
        F_unif=vector(mode = "double", length = reference_size)
        f_unif=vector(mode = "double", length = reference_size)
        F_log=vector(mode = "double", length = reference_size)
        f_log=vector(mode = "double", length = reference_size)
        
        grids_full=.C("get_grids_full",
        file=full_data,
        columns=as.integer(columns_full_data),
        full_tests=as.double(full_tests),
        rows=as.double(full_tests),
        stat_column=as.integer(full_column_id),
        grid_size=as.integer(reference_size),
        grid_unif=grid_unif,
        grid_log=grid_log,
        F_log=F_log,
        f_log=f_log,
        F_unif=F_unif,
        f_unif=f_unif,
        min_stat=min_stat,
        max_stat=max_stat,
        zero_count=as.integer(0),
        print_messages=print_messages,
        rownames_flag=as.integer(rownames_flag_full))
        
        grid_unif<-grids_full$grid_unif
        grid_log<-grids_full$grid_log
        f_log<-grids_full$f_log
        F_log<-grids_full$F_log
        f_unif<-grids_full$f_unif
        F_unif<-grids_full$F_unif
        zero_count_full<-grids_full$zero_count
        rm(grids_full)
        
        if(flag_debug_august==TRUE){
            df=data.frame(grid_unif,f_unif,F_unif,grid_log,f_log,F_log)
            write.table(df,paste0(output_name,"/file_full_dist.txt"))
        }
        
        log_max_ratio=log10(max_stat/min_stat)
        delta_log=log_max_ratio/(reference_size-1)
        
    }else{
        
        full_data[,full_column_id]=as.double(as.matrix(full_data[,full_column_id]))
        if(length(full_data)>1){full_data=full_data[order(full_data[,full_column_id]),]}
        else{full_data=as.matrix(full_data[order(full_data[,full_column_id]),])}
        zero_count_full=length(which(full_data[,full_column_id]==0))
        data=full_data[1:full_tests,full_column_id]
        grid=unique(data)
        
        reference_size=length(grid)
        
        F=vector(mode = "double", length = reference_size)
        f=vector(mode = "double", length = reference_size)
        min_stat=grid[1]
        max_stat=grid[reference_size]
        
        dist=.C("get_distribution_df",
        data=as.double(data),
        data_size=as.integer(full_tests),
        processed_data_size=as.integer(full_tests),
        grid=as.double(grid),
        grid_size=as.integer(reference_size),
        f=f,
        F=F)
        
        f<-dist$f
        F<-dist$F
        rm(dist)
        
        if(flag_debug_august==TRUE){
            df=data.frame(grid,f,F)
            write.table(df,paste0(output_name,"/df_full_dist.txt"))
            
        }
    }
    
    if(print_messages)
    {
        message("Obtaining p-values distribution from permutation tests...")
    }
    
    if(perm_data_is_file)
    {
        if(full_data_is_file)
        {
            F_o_unif=vector(mode = "double", length = reference_size)
            f_o_unif=vector(mode = "double", length = reference_size)
            F_o_log=vector(mode = "double", length = reference_size)
            f_o_log=vector(mode = "double", length = reference_size)
            
            grids_perm=.C("get_grids_perm",
            file=perm_data,
            columns=as.integer(columns_perm_data),
            perm_tests=as.double(perm_tests),
            rows=as.double(perm_tests),
            stat_columns=length(perm_column_ids),
            stat_columns_array=as.integer(perm_column_ids),
            grid_size=as.integer(reference_size),
            grid_unif=grid_unif,
            grid_log=grid_log,
            F_o_log=F_o_log,
            f_o_log=f_o_log,
            F_o_unif=F_o_unif,
            f_o_unif=f_o_unif,
            zero_count=as.integer(0),
            print_messages=print_messages,
            rownames_flag=as.integer(rownames_flag_perm))
            
            f_o_log<-grids_perm$f_o_log
            F_o_log<-grids_perm$F_o_log
            f_o_unif<-grids_perm$f_o_unif
            F_o_unif<-grids_perm$F_o_unif
            zero_count_perm<-grids_perm$zero_count
            rm(grids_perm)
            
            if(flag_debug_august==TRUE){
                df=data.frame(grid_unif,f_o_unif,F_o_unif,grid_log,f_o_log,F_o_log)
                write.table(df,paste0(output_name,"/file_file_perm_dist.txt"))
            }
        }else{
            F_o=vector(mode = "double", length = reference_size)
            f_o=vector(mode = "double", length = reference_size)
            
            grids_perm=.C("get_grids_perm_from_data",
            file=perm_data,
            columns=as.integer(columns_perm_data),
            perm_tests=as.double(perm_tests),
            rows=as.double(perm_tests),
            stat_columns=length(perm_column_ids),
            stat_columns_array=as.integer(perm_column_ids),
            grid_size=as.integer(reference_size),
            grid=grid,
            F_o=F_o,
            f_o=f_o,
            zero_count=as.integer(0),
            print_messages=print_messages,
            rownames_flag=as.integer(rownames_flag_perm))
            
            f_o<-grids_perm$f_o
            F_o<-grids_perm$F_o
            zero_count_perm<-grids_perm$zero_count
            rm(grids_perm)
            
            if(flag_debug_august==TRUE){
                df=data.frame(grid,f_o,F_o)
                write.table(df,paste0(output_name,"/df_file_perm_dist.txt"))
            }
        }
    }else{
        perm_data=perm_data[order(perm_data)]
        
        if(full_data_is_file)
        {
            F_o_log=vector(mode = "double", length = reference_size)
            f_o_log=vector(mode = "double", length = reference_size)
            F_o_unif=vector(mode = "double", length = reference_size)
            f_o_unif=vector(mode = "double", length = reference_size)
            
            dist_perm=.C("get_grids_perm_df",
            perm_tests_array=as.double(as.matrix(perm_data)),
            perm_tests=as.integer(perm_tests),
            perm_tests=as.integer(perm_tests),
            grid_size=as.integer(reference_size),
            grid_log=grid_log,
            grid_unif=grid_unif,
            F_o_log=F_o_log,
            f_o_log=f_o_log,
            F_o_unif=F_o_unif,
            f_o_unif=f_o_unif,
            zero_count=as.integer(0))
            
            F_o_log=as.double(dist_perm$F_o_log)
            F_o_unif=as.double(dist_perm$F_o_unif)
            f_o_log=as.double(dist_perm$f_o_log)
            f_o_unif=as.double(dist_perm$f_o_unif)
            zero_count_perm<-dist_perm$zero_count
            rm(dist_perm)
            
            if(flag_debug_august==TRUE){
                df=data.frame(grid_unif,f_o_unif,F_o_unif,grid_log,f_o_log,F_o_log)
                write.table(df,paste0(output_name,"/file_df_perm_dist.txt"))
            }
            
        }else{
            
            zero_count_perm=length(which(perm_data==0))
            F_o=vector(mode = "double", length = reference_size)
            f_o=vector(mode = "double", length = reference_size)
            
            dist_perm=.C("get_distribution_df",
            data=as.double(perm_data),
            data_size=as.integer(perm_tests),
            processed_data_size=as.integer(perm_tests),
            grid=as.double(grid),
            grid_size=as.integer(reference_size),
            f=f_o,
            F=F_o)
            
            f_o<-dist_perm$f
            F_o<-dist_perm$F
            rm(dist_perm)
            
            if(flag_debug_august==TRUE){
                df=data.frame(grid,f_o,F_o)
                write.table(df,paste0(output_name,"/df_df_perm_dist.txt"))
            }
        }
    }
    
    if(print_messages)
    {
        message("\n_____________________ STEP 3: Obtain fraction of true negatives pi_o ________________________________________")
    }
    
    F_g=vector(mode = "double", length = reference_size)
    F_o_g=vector(mode = "double", length = reference_size)
    f_g=vector(mode = "double", length = reference_size)
    f_o_g=vector(mode = "double", length = reference_size)
    
    fitting_goodness=vector(mode = "double", length = 2)
    flag_grenander_choice=vector(mode = "integer", length = 2)
    
    if(full_data_is_file){
        grid=grid_unif
        F=F_unif
        F_o=F_o_unif
        f=f_unif
        f_o=f_o_unif
    }

    rho_o_estimation_list=estimate_rho_o(grid,F,f,F_o,f_o,print_messages,reference_size,p_star_exhaustive_sweep,
    full_data_is_file,rho_o_splines_knots,rho_o_splines_knots_range,rho_o_splines_iterations)
    if(print_messages)
    {
    message("Obtaining pi_o...")
    }
    index_p_star=rho_o_estimation_list$index_p_star
    p_star=rho_o_estimation_list$p_star
    rho_o_scaling=rho_o_estimation_list$rho_o_scaling
    flag_grenander_choice=rho_o_estimation_list$flag_grenander_choice
    fitting_goodness=rho_o_estimation_list$fitting_goodness
    
    if(print_messages)
    {
        message("\n_____________________ STEP 4: Obtain false discovery rates: permFDR method __________________________________")
    }
    
    if(full_data_is_file)
    {
        index_p_star_unif=index_p_star
        index_p_star_log=(as.integer)(log10(p_star/min_stat)/delta_log+1)
        
        grid=grid_log
        F=F_log
        F_o=F_o_log
        index_p_star=index_p_star_log
    }
    
    pi_o=1+rho_o_scaling*F_o[index_p_star]-F[index_p_star]
    pi_o=max(min(pi_o,1),0)
    
    F_null=F_o*rho_o_scaling
    F_null[index_p_star:reference_size]=F_o[index_p_star]*rho_o_scaling+F[index_p_star:reference_size]-F[index_p_star]
    
    f_null=f_o*rho_o_scaling
    f_null[index_p_star:reference_size]=f[index_p_star:reference_size]

    Fdr=F_null/F
    Fndr=(1-F-pi_o+F_null)/(1-F)
    Fndr[length(Fndr)]=Fndr[length(Fndr)-1]
    
    zero_Fdr=rep(1,8)
    if(zero_count_full!=0 & zero_count_full!=full_tests)
    {
        zero_Fdr[1]=pi_o*zero_count_perm/zero_count_full
        zero_Fdr[2]=((full_tests-zero_count_full)/full_tests -pi_o +rho_o_scaling*zero_count_perm/perm_tests)*full_tests/(full_tests-zero_count_full)
    }else if(zero_count_full==0){
        zero_Fdr[1]=pi_o
        zero_Fdr[2]=min(1,1 -pi_o +rho_o_scaling*zero_count_perm/perm_tests)
    }else{
        zero_Fdr[1]=pi_o*zero_count_perm/zero_count_full
        zero_Fdr[2]=1
    }
        
    alt_methods_aux=as.integer(valid_alt_methods %in% alt_methods)
    
    #### Comprobar que f_null<f###
    
    maximum_scaling_error_bubble=check_scaling(grid,f_null,f)
    if(maximum_scaling_error_bubble$scaling_error>0.05){
    message("Warning: f_null overpasses f by an integral fraction (>5%) equal to: ",maximum_scaling_error_bubble$scaling_error," between p=",maximum_scaling_error_bubble$end," and p=",maximum_scaling_error_bubble$end)
    }
    ###################
    
    Fdr_ST=0
    Fndr_ST=0
    Fdr_BH_unif=0
    Fdr_BH_perm=0
    Fdr_MV=0
    BF=0
    
    if(length(which(alt_methods_aux==1))>0 & print_messages)
    {
        message("\n_____________________ STEP 5: Obtain false discovery rates: Alternative methods _____________________________")
    }
    
    if(full_data_is_file)
    {
		F=F_log
		F_o=F_o_log
		grid=grid_log
    }
        
    if(length(which(alt_methods_aux==1))>0)
    {
        maximum_scaling_error_bubble=c(0,0,0)
        if("Fdr_ST" %in% alt_methods)
        {
            if(print_messages){message("Calculating Storey-Tibshinari FDRs...")}
            
				if(full_data_is_file)
				{
					grid=grid_unif
					F=F_unif
					F_o=F_o_unif
				}
                
				rho_hat_ST=(1-F)/(1-F_o)
				finite_set_ST=which(is.finite(as.numeric(as.character(rho_hat_ST))))
                
            
                constraint_matrix=as.matrix(data.frame(c(0,2),c(0,grid[finite_set_ST[length(finite_set_ST)]]),c(1,0)))


				f_hat_ST=splines_fit(grid[finite_set_ST],rho_hat_ST[finite_set_ST],rho_o_splines_knots,rho_o_splines_iterations,
                rho_o_splines_knots_range,trend="decrease",full_number=1,perm_number=1,print_messages=print_messages,constraint_matrix=constraint_matrix)
            
				pi_hat_fitted_ST<-f_hat_ST$fitted
				pi_o_ST=pi_hat_fitted_ST[length(pi_hat_fitted_ST)]
				pi_o_ST=max(min(pi_o_ST,1),0)
				rm(f_hat_ST)
                
				if(full_data_is_file)
				{
					grid=grid_log
					F=F_log
					F_o=F_o_log
				}
                f_null=rep(pi_o_ST,length(f))
                maximum_scaling_error_bubble=check_scaling(grid,f_null,f)
                if(maximum_scaling_error_bubble$scaling_error>0.05){
                    message("Warning: Storey-Tibshirani method: f_null overpasses f by an integral fraction (>5%) equal to: ",maximum_scaling_error_bubble$scaling_error," between p=",maximum_scaling_error_bubble$end," and p=",maximum_scaling_error_bubble$end)
                }
                Fdr_ST=pi_o_ST*grid/F
                if(zero_count_full!=0)
                {
                    zero_Fdr[3]=pi_o_ST*zero_count_perm/zero_count_full
                }else{
                    zero_Fdr[3]=pi_o
                }
        }
        if("Fndr_ST" %in% alt_methods)
        {
            if(print_messages){message("Calculating Storey-Tibshinari FNDRs...")}
            
            if(!("Fdr_ST" %in% alt_methods))
            {
                if(full_data_is_file)
                {
                    grid=grid_unif
                    F=F_unif
                    F_o=F_o_unif
                }
                
                rho_hat_ST=(1-F)/(1-F_o)
                finite_set_ST=which(is.finite(as.numeric(as.character(rho_hat_ST))))
                constraint_matrix=as.matrix(data.frame(c(0,2),c(0,grid[finite_set_ST[length(finite_set_ST)]]),c(1,0)))
            f_hat_ST=splines_fit(grid[finite_set_ST],rho_hat_ST[finite_set_ST],rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,trend="decrease",full_number=1,perm_number=1,print_messages=print_messages,constraint_matrix=constraint_matrix)
                
                rho_hat_fitted_ST<-f_hat_ST$fitted
                if(flag_debug_august==TRUE){
                    df=data.frame(x=grid[finite_set_ST],y=rho_hat_ST[finite_set_ST],yfit=rho_hat_fitted_ST)
                    p=ggplot(df)+geom_point(aes(x=x,y=y))+geom_line(aes(x=x,y=yfit))
                    pdf("ST_when_entering_by_fndr_fit.pdf")
                    print(p)
                    dev.off()
                }
                pi_o_ST=rho_hat_fitted_ST[length(rho_hat_fitted_ST)]
                pi_o_ST=min(pi_o_ST,1)
                pi_o_ST=max(pi_o_ST,0)
                rm(f_hat_ST)
                
                if(full_data_is_file)
                {
                    grid=grid_log
                    F=F_log
                    F_o=F_o_log
                }
                f_null=rep(pi_o_ST,length(f))
                maximum_scaling_error_bubble=c(0,0,0)
                maximum_scaling_error_bubble=check_scaling(grid,f_null,f)
                if(maximum_scaling_error_bubble$scaling_error>0.05){
                    message("Warning: Storey-Tibshirani method: f_null overpasses f by an integral fraction (>5%) equal to: ",maximum_scaling_error_bubble$scaling_error," between p=",maximum_scaling_error_bubble$end," and p=",maximum_scaling_error_bubble$end)
                }
            }
            Fndr_ST=(1-F-pi_o_ST*(1-grid))/(1-F)
            Fndr_ST[length(Fndr_ST)]=Fndr_ST[length(Fndr_ST)-1]
            if(zero_count_full!=full_tests)
            {
                zero_Fdr[4]=((full_tests-zero_count_full)/full_tests -pi_o_ST*(1-zero_count_perm/perm_tests))*full_tests/(full_tests-zero_count_full)
            }else{
                zero_Fdr[4]=1
            }
		}
        if("Fdr_BH_unif" %in% alt_methods)
        {
            if(print_messages){message("Calculating Benjamini-Hochberg FDRs (uniform null)...")}
            Fdr_BH_unif=grid/F
            zero_Fdr[5]=0
            
        }
        if("Fdr_BH_perm" %in% alt_methods)
        {
            if(print_messages){message("Calculating Benjamini-Hochberg FDRs (permutations-based null)...")}
            Fdr_BH_perm=F_o/F
            if(zero_count_full!=0)
            {
                zero_Fdr[6]=zero_count_perm/zero_count_full
            }else{
                zero_Fdr[6]=1
            }
        }
        if("Fdr_MV" %in% alt_methods)
        {
            maximum_scaling_error_bubble=c(0,0,0)
            if(print_messages){message("Calculating Millstein-Volfson FDRs...")}
            rho_hat_MV=(1-F)/(1-F_o)
            end_plateau_MV=which(F_o==1)
            rho_hat_MV[end_plateau_MV]=rho_hat_MV[end_plateau_MV[1]-1]
            
            ## In MV, an adequate fit along all the range of p-values is more important than in permFDR, where only the final trend of the fit defines rho_o.
            ## For that reason, if the signal is extremely strong, an alternative fitting strategy is implemented, with two distinctive features:
            ## 1) Increased number of knots.
            ## 2) x-coordinates are log transform to allow a better resolution of the highly popolated
            
            grid_aux=grid[which(grid<1e-50)]
            if(length(grid_aux)>0){
                peak=F[length(grid_aux)]
            }else{peak=0}
            if(peak>0.10) #Si hay un 10% de los datos por debajo de 1e-50
            {
                rho_o_splines_knots=min(30,as.integer(2+9.3*log(length(grid)/10)))
                if(grid[1]==0){
                    grid_aux=grid[2:length(grid)]
                    rho_hat_MV_aux=rho_hat_MV[2:length(grid)]
                    constraint_matrix=as.matrix(data.frame(c(0),c(log10(grid_aux)[1]),c(rho_hat_MV_aux[1])))
                    f_hat_MV_aux=splines_fit(log10(grid_aux),rho_hat_MV_aux,rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,trend="none",full_number=1,perm_number=1,print_messages=print_messages,constraint_matrix=constraint_matrix)
                    rho_hat_fitted_MV_aux<-f_hat_MV_aux$fitted
                    rho_hat_fitted_MV=c(rho_hat_MV[1],rho_hat_fitted_MV_aux)
                    rho_hat_long<-rho_hat_fitted_MV
                }else{
                    constraint_matrix=as.matrix(data.frame(c(0),c(log10(grid)[1]),c(rho_hat_MV[1])))
                    f_hat_MV=splines_fit(log10(grid),rho_hat_MV,rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,trend="none",full_number=1,perm_number=1,print_messages=print_messages,constraint_matrix=constraint_matrix)
                    rho_hat_fitted_MV<-f_hat_MV$fitted
                    rho_hat_long<-rho_hat_fitted_MV
                }
            }else{
                f_hat_MV=splines_fit(grid,rho_hat_MV,rho_o_splines_knots,rho_o_splines_iterations,rho_o_splines_knots_range,trend="none",full_number=1,perm_number=1,print_messages=print_messages)
                rho_hat_fitted_MV<-f_hat_MV$fitted
                rho_hat_long<-rho_hat_fitted_MV
            }
            
            if(flag_debug_august==TRUE){
                df=data.frame(x=grid,y=rho_hat_MV,yfit=rho_hat_long)
                write.table(df,paste0(output_name,"/MV_fit_output.txt"))
                p=ggplot(df)+geom_point(aes(x=x,y=y))+geom_line(aes(x=x,y=yfit))
                pp=ggplot(df)+geom_point(aes(x=log10(x),y=y))+geom_line(aes(x=log10(x),y=yfit))
                
                pdf("MV_fit.pdf")
                print(p)
                dev.off()
                pdf("MV_fit_log.pdf")
                print(pp)
                dev.off()
            }
            f_null=rho_hat_long*f_o
            maximum_scaling_error_bubble=check_scaling(grid,f_null,f)
            if(maximum_scaling_error_bubble$scaling_error>0.05){
                message("Warning: Millstein-Volfson method: f_null overpasses f by an integral fraction (>5%) equal to: ",maximum_scaling_error_bubble$scaling_error," between p=",maximum_scaling_error_bubble$end," and p=",maximum_scaling_error_bubble$end)
            }
            Fdr_MV=rho_hat_long*F_o/F
            if(zero_count_full!=0)
            {
                zero_Fdr[7]=rho_hat_long[1]*zero_count_perm/zero_count_full
            }else{
                zero_Fdr[7]=rho_hat_long[1]
            }
        }
        if("BF" %in% alt_methods)
        {
            if(print_messages){message("Calculating Bonferroni adjusted p-values...")}
            
            BF=grid*length(grid)
            BF[which(BF>1)]=1
            zero_Fdr[8]=0
        }
    }
    
    if(full_data_is_file){
        significant_tests=vector(mode = "double", length = (2+sum(alt_methods_aux)))
        if(print_messages){message("Interpolating FDRs from grid p-values...")}
        result=.C("interpolate_fdrs",
        file=full_data,
        columns=as.integer(columns_full_data),
        flag_rownames=as.integer(rownames_flag_full),
        rows=as.double(full_tests),
        stat_column=as.integer(full_column_id),
        sampled_stats=as.integer(reference_size),
        sampled_stats_array=grid_log,
        Fdr=Fdr,
        Fndr=Fndr,
        Fdr_ST=Fdr_ST,
        Fndr_ST=Fndr_ST,
        Fdr_BH_unif=Fdr_BH_unif,
        Fdr_BH_perm=Fdr_BH_perm,
        Fdr_MV=Fdr_MV,
        BF=BF,
        name=output_name,
        alt_methods=as.integer(alt_methods_aux),
        zero_Fdr=zero_Fdr,
        report_threshold=as.double(report_threshold),
        print_messages=print_messages,
        significance_threshold=significance_threshold,
        significant_tests=as.double(significant_tests))
        significant_tests_vector=as.double(result$significant_tests)
        
        Fdr<-as.double(result$Fdr)
        Fndr<-as.double(result$Fndr)
        Fdr_ST<-as.double(result$Fdr_ST)
        Fndr_ST<-as.double(result$Fndr_ST)
        Fdr_BH_unif<-as.double(result$Fdr_BH_unif)
        Fdr_BH_perm<-as.double(result$Fdr_BH_perm)
        Fdr_MV<-as.double(result$Fdr_MV)
        BF<-as.double(result$BF)
        
        significant_tests=data.frame(significant_tests_vector[1],significant_tests_vector[2])
        cont=1
        for(i in 1:length(alt_methods_aux))
        {
            if(alt_methods_aux[i]==1)
            {
                significant_tests=cbind(significant_tests,significant_tests_vector[cont+2])
                cont=cont+1
            }
        }
        if(sum(alt_methods_aux)>0){colnames(significant_tests)=c("Fdr_SAB","Fndr_SAB",alt_methods)}
        else{colnames(significant_tests)=c("Fdr_SAB","Fndr_SAB")}
    }else{
        result=.C("monotone_fdrs",Fdr=Fdr,Fndr=Fndr,
        Fdr_ST=Fdr_ST,
        Fndr_ST=Fndr_ST,
        Fdr_BH_unif=Fdr_BH_unif,
        Fdr_BH_perm=Fdr_BH_perm,
        Fdr_MV=Fdr_MV,
        alt_methods=as.integer(alt_methods_aux),
        sampled_stats=as.integer(reference_size))
        
        Fdr<-as.double(result$Fdr)
        Fndr<-as.double(result$Fndr)
        Fdr_ST<-as.double(result$Fdr_ST)
        Fndr_ST<-as.double(result$Fndr_ST)
        Fdr_BH_unif<-as.double(result$Fdr_BH_unif)
        Fdr_BH_perm<-as.double(result$Fdr_BH_perm)
        Fdr_MV<-as.double(result$Fdr_MV)
        
        significant_tests=data.frame(0,0)
        if(sum(alt_methods_aux)>0){
            for(i in 1:sum(alt_methods_aux))
            {
                significant_tests=cbind(significant_tests,0)
            }
        }
        if(sum(alt_methods_aux)>0){colnames(significant_tests)=c("Fdr_SAB","Fndr_SAB",alt_methods)}
        else{colnames(significant_tests)=c("Fdr_SAB","Fndr_SAB")}
        
        significant_tests$Fdr_SAB=length(which(Fdr<=significance_threshold))
        significant_tests$Fndr_SAB=length(which(Fndr<=significance_threshold))
        if("Fdr_ST" %in% alt_methods){  significant_tests$Fdr_ST=length(which(Fdr_ST<=significance_threshold))  }
        if("Fndr_ST" %in% alt_methods){  significant_tests$Fndr_ST=length(which(Fndr_ST<=significance_threshold))  }
        if("Fdr_BH_unif" %in% alt_methods){  significant_tests$Fdr_BH_unif=length(which(Fdr_BH_unif<=significance_threshold))  }
        if("Fdr_BH_perm" %in% alt_methods){  significant_tests$Fdr_BH_perm=length(which(Fdr_BH_perm<=significance_threshold))  }
        if("Fdr_MV" %in% alt_methods){  significant_tests$Fdr_MV=length(which(Fdr_MV<=significance_threshold))  }
        if("BF" %in% alt_methods){  significant_tests$BF=length(which(BF<=significance_threshold))  }
    }
    
    if(print_messages)
    {
        if(length(which(alt_methods_aux==1))>0)
        {
            message("\n_____________________ STEP 6: Print outputs _________________________________________________________________")
        }else{
            message("\n_____________________ STEP 5: Print outputs _________________________________________________________________")
        }
    }
    
    output=list(full_tests,perm_tests,rho_o_scaling,pi_o,p_star,significant_tests,call=deparse(match.call()))
    names(output)=c("full_tests","perm_tests","rho_o","pi_o","p_star","hits","call")
    
    if("Fdr_ST" %in% alt_methods | "Fndr_ST" %in% alt_methods)
    {
        names=c(names(output),"pi_o_ST")
        output[[length(output)+1]]=pi_o_ST
        names(output)=names
    }
    
    if(full_data_is_file==0)
    {
        fdrs=data.frame(full_data,Fdr_SAB=-1,Fndr_SAB=-1)
        first_entries=which(!duplicated(fdrs[,full_column_id]))
        duplicated_entries=which(duplicated(fdrs[,full_column_id]))
        
        fdrs$Fdr_SAB[first_entries]=Fdr
        fdrs$Fndr_SAB[first_entries]=Fndr
        
        for(i in 1:length(duplicated_entries))
        {
            fdrs$Fdr_SAB[duplicated_entries[i]]=fdrs$Fdr_SAB[duplicated_entries[i]-1]
            fdrs$Fndr_SAB[duplicated_entries[i]]=fdrs$Fndr_SAB[duplicated_entries[i]-1]
        }
        
        if(alt_methods_aux[1]==1)
        {
            fdrs=data.frame(fdrs,Fdr_ST=0)
            fdrs$Fdr_ST[first_entries]=Fdr_ST
            
            for(i in 1:length(duplicated_entries))
            {
                fdrs$Fdr_ST[duplicated_entries[i]]=fdrs$Fdr_ST[duplicated_entries[i]-1]
            }
        }
        
        if(alt_methods_aux[2]==1)
        {
            fdrs=data.frame(fdrs,Fndr_ST=0)
            fdrs$Fndr_ST[first_entries]=Fndr_ST
            for(i in 1:length(duplicated_entries))
            {
                fdrs$Fndr_ST[duplicated_entries[i]]=fdrs$Fndr_ST[duplicated_entries[i]-1]
            }
        }
        
        if(alt_methods_aux[3]==1)
        {
			fdrs=data.frame(fdrs,Fdr_BH_unif=0)
            fdrs$Fdr_BH_unif[first_entries]=Fdr_BH_unif
            for(i in 1:length(duplicated_entries))
            {
                fdrs$Fdr_BH_unif[duplicated_entries[i]]=fdrs$Fdr_BH_unif[duplicated_entries[i]-1]
            }
        }
        
        if(alt_methods_aux[4]==1)
        {
			fdrs=data.frame(fdrs,Fdr_BH_perm=0)
            fdrs$Fdr_BH_perm[first_entries]=Fdr_BH_perm
            for(i in 1:length(duplicated_entries))
            {
                fdrs$Fdr_BH_perm[duplicated_entries[i]]=fdrs$Fdr_BH_perm[duplicated_entries[i]-1]
            }
        }
        
        if(alt_methods_aux[5]==1)
        {
            fdrs=data.frame(fdrs,Fdr_MV=0)
            fdrs$Fdr_MV[first_entries]=Fdr_MV
            for(i in 1:length(duplicated_entries))
            {
                fdrs$Fdr_MV[duplicated_entries[i]]=fdrs$Fdr_MV[duplicated_entries[i]-1]
            }
        }
        
        if(alt_methods_aux[6]==1)
        {
            fdrs=data.frame(fdrs,BF=0)
            fdrs$BF[first_entries]=BF
            for(i in 1:length(duplicated_entries))
            {
                fdrs$BF[duplicated_entries[i]]=fdrs$BF[duplicated_entries[i]-1]
            }
        }
        
        names=c(names(output),"fdrs")
        output[[length(output)+1]]=fdrs
        names(output)=names
    }
    
   if(plot==TRUE)
    {
        ###Build data.frames to plot: 1) local distributions 2) Grenander (local) distributions. 3) ECDFs 4)fdrs
        
        f_g=rho_o_estimation_list$f_g
        F_g=rho_o_estimation_list$F_g
        f_o_g=rho_o_estimation_list$f_o_g
        F_o_g=rho_o_estimation_list$F_o_g
        rho_fit_table=rho_o_estimation_list$rho_fit_table
        
        if(p_star_exhaustive_sweep==TRUE)
        {
            rho_o_scaling_loop_table=rho_o_estimation_list$rho_o_scaling_loop_table
        }

        if(full_data_is_file)
        {
            grid=grid_log
            F=F_log
        }
        
        
        ##1. Load fdrs table
        fdrs=data.frame(statistic=grid,quantile=(F/F[length(F)])*full_tests,Fdr=Fdr,Fndr=Fndr)
        labels_fdr=c(expression(Fdr[{SAB}]),expression(Fndr[{SAB}]))
        
        if(alt_methods_aux[1]==1)
        {
            fdrs=data.frame(fdrs,Fdr_ST=Fdr_ST)
            labels_fdr=c(labels_fdr,expression(Fdr[{ST}]))
        }
        
        if(alt_methods_aux[2]==1)
        {
            fdrs=data.frame(fdrs,Fndr_ST=Fndr_ST)
            labels_fdr=c(labels_fdr,expression(Fndr[{ST}]))
        }
        
        if(alt_methods_aux[3]==1)
        {
            fdrs=data.frame(fdrs,Fdr_BH_unif=Fdr_BH_unif)
            labels_fdr=c(labels_fdr,expression(Fdr[{BH(unif)}]))
        }
        
        if(alt_methods_aux[4]==1)
        {
            fdrs=data.frame(fdrs,Fdr_BH_perm=Fdr_BH_perm)
            labels_fdr=c(labels_fdr,expression(Fdr[{BH(perm)}]))
        }
        
        if(alt_methods_aux[5]==1)
        {
            fdrs=data.frame(fdrs,Fdr_MV=Fdr_MV)
            labels_fdr=c(labels_fdr,expression(Fdr[{MV}]))
        }
        
        if(alt_methods_aux[6]==1)
        {
            fdrs=data.frame(fdrs,BF=BF)
            labels_fdr=c(labels_fdr,expression(BF))
        }
        
        ## Restrict it to the segment of interest
        if(report_threshold<max_stat)
        {
            end=fdrs$statistic[which(fdrs$Fdr>=report_threshold)][1]
            if("Fdr_ST" %in% alt_methods) end=max(end,fdrs$statistic[which(fdrs$Fdr_ST>=report_threshold)][1])
            if("Fdr_BH_perm" %in% alt_methods) end=max(end,fdrs$statistic[which(fdrs$Fdr_BH_perm>=report_threshold)][1])
            if("Fdr_BH_unif" %in% alt_methods) end=max(end,fdrs$statistic[which(fdrs$Fdr_BH_unif>=report_threshold)][1])
            if("Fdr_MV" %in% alt_methods) end=max(end,fdrs$statistic[which(fdrs$Fdr_MV>=report_threshold)][1])
            if(is.na(end)) end=max_stat
        }else{
            end=max_stat
        }
        
        ##Load tab: dataframe with distributions, both local and cumulative (some are still provisional)
        if(full_data_is_file==0)
        {
            F_null_g=F_o_g*rho_o_scaling
            F_null_g[index_p_star:reference_size]=(F_o_g[index_p_star]*rho_o_scaling+F_g[index_p_star:reference_size]-F_g[index_p_star])
            pi_o_g=F_null_g[length(F_null_g)]
            F_null_g=F_null_g/F_null_g[length(F_null_g)]
            tab<-data.frame(statistic=grid,f_o=f_o,f=f,f_o_g=f_o_g,f_g=f_g,F_o=F_o,F=F,F_o_g=F_o_g,F_g=F_g,pi_o_F_null=pi_o_g*F_null_g,pi_o_f_null=pi_o_g*F_null_g)
        }else{
            F_null_unif=F_o_g*rho_o_scaling
            F_null_unif[index_p_star_unif:reference_size]=(F_o_g[index_p_star_unif]*rho_o_scaling+F_g[index_p_star_unif:reference_size]-F_g[index_p_star_unif])
            pi_o_g=F_null_unif[length(F_null_unif)]
            F_null_unif=F_null_unif/F_null_unif[length(F_null_unif)]
            tab<-data.frame(statistic=grid_unif,f_o=f_o_unif,f=f_unif,f_o_g=f_o_g,f_g=f_g,F_o=F_o_unif,F=F_unif,F_o_g=F_o_g,F_g=F_g,pi_o_F_null=pi_o_g*F_null_unif,pi_o_f_null=pi_o_g*F_null_unif)
        }
        
        cut=which(tab$statistic>=end*0.999999)[1]
        tab=tab[1:cut,]
        
        
        ###Reduce points/bins to plot to a constant number (to avoid generating pdfs with as many points as tests, which might be a lot) and get correct versions of the distributions. Report threshold will be accepted only if allows plotting at least 20 points.
        
        bins_array=c(500,450,400,350,300,250,200,150,100,50,20)
        enough_bins=0;
        for(bins_number in c(500,450,400,350,300,250,200,150,100,50,20))
        {
            if(nrow(tab)>bins_number)
            {
                enough_bins=1
                break
            }else{
                bins_array=bins_array[(which(bins_array==bins_number)+1):length(bins_array)]
            }
        }
        if(enough_bins==0)
        {
            end=max_stat
            if(full_data_is_file){
                tab<-data.frame(statistic=grid_unif,f_o=f_o_unif,f=f_unif,f_o_g=f_o_g,f_g=f_g,F_o=F_o_unif,F=F_unif,F_o_g=F_o_g,F_g=F_g,pi_o_F_null=pi_o_g*F_null_unif,pi_o_f_null=pi_o_g*F_null_unif)
            }else{
                tab<-data.frame(statistic=grid,f_o=f_o,f=f,f_o_g=f_o_g,f_g=f_g,F_o=F_o,F=F,F_o_g=F_o_g,F_g=F_g,pi_o_F_null=pi_o_g*F_null_g,pi_o_f_null=pi_o_g*F_null_g)
            }
            cut=which(tab$statistic>=end*0.999999)[1]
            tab=tab[1:cut,]
            message("WARNING: Not enough resolution for this report_threshold. Set to 1 in the plot")
        }
        
        for(bins_number in bins_array)
        {
            if(full_data_is_file)
            {
                reduction_factor=as.integer(nrow(tab)/bins_number)
                v_aux=seq(1,nrow(tab),reduction_factor)
                v_aux=c(v_aux[2:length(v_aux)],nrow(tab))
                v_aux=unique(v_aux)
            }else{
                reduction_factor=as.integer(nrow(tab)/bins_number)
                v_aux=vector(mode = "double", length = bins_number)
                delta=as.double((end-min_stat)/bins_number)
                grid_aux=seq(min_stat,end,delta)
                
                index_aux=1
                for(i in 2:length(v_aux))
                {
                    while(grid_aux[i]>grid[index_aux])
                    index_aux=index_aux+1
                    v_aux[i]=index_aux
                }
                v_aux=c(v_aux[2:length(v_aux)],nrow(tab))
                v_aux=unique(v_aux)
            }
            
            v=v_aux
            rm(v_aux)
            
            tab_reduced=tab[v,]
            tab_reduced$f_o[1]=tab_reduced$F_o[1]/(tab_reduced$statistic[1]-min_stat)
            tab_reduced$f_o_g[1]=tab_reduced$F_o_g[1]/(tab_reduced$statistic[1]-min_stat)
            tab_reduced$f[1]=tab_reduced$F[1]/(tab_reduced$statistic[1]-min_stat)
            tab_reduced$f_g[1]=tab_reduced$F_g[1]/(tab_reduced$statistic[1]-min_stat)
            tab_reduced$pi_o_f_null[1]=tab_reduced$pi_o_F_null[1]/(tab_reduced$statistic[1]-min_stat)
            
            for(i in 2:nrow(tab_reduced))
            {
                tab_reduced$f_o[i]=(tab_reduced$F_o[i]-tab_reduced$F_o[i-1])/(tab_reduced$statistic[i]-tab_reduced$statistic[i-1])
                tab_reduced$f_o_g[i]=(tab_reduced$F_o_g[i]-tab_reduced$F_o_g[i-1])/(tab_reduced$statistic[i]-tab_reduced$statistic[i-1])
                tab_reduced$f[i]=(tab_reduced$F[i]-tab_reduced$F[i-1])/(tab_reduced$statistic[i]-tab_reduced$statistic[i-1])
                tab_reduced$f_g[i]=(tab_reduced$F_g[i]-tab_reduced$F_g[i-1])/(tab_reduced$statistic[i]-tab_reduced$statistic[i-1])
                tab_reduced$pi_o_f_null[i]=(tab_reduced$pi_o_F_null[i]-tab_reduced$pi_o_F_null[i-1])/(tab_reduced$statistic[i]-tab_reduced$statistic[i-1])
            }
            maxim=max(tab_reduced$f_o,tab_reduced$f,tab_reduced$f_o_g,tab_reduced$f_g,tab_reduced$pi_o_f_null)
            if(maxim<40)
            {break}
        }
        tab_reduced=tab_reduced[-nrow(tab_reduced),]
        tab_reduced=rbind(tab_reduced[1,],tab_reduced)
        tab_reduced[1,1]=min_stat
        
        tab_1=data.frame(tab_reduced[,c(1,2)],Data="Permuted")
        colnames(tab_1)=c("statistic","series","Data")
        
        tab_2=data.frame(tab_reduced[,c(1,3)],Data="Full")
        colnames(tab_2)=c("statistic","series","Data")
        
        tab_hists=rbind(tab_1,tab_2)
        
        tab_3=data.frame(tab_reduced[,c(1,4)],Estimated_distribution="f_perm")
        colnames(tab_3)=c("statistic","series","Estimated_distribution")
        
        tab_4=data.frame(tab_reduced[,c(1,5)],Estimated_distribution="f_full")
        colnames(tab_4)=c("statistic","series","Estimated_distribution")
        
        tab_5=data.frame(tab_reduced[,c(1,11)],Estimated_distribution="pi_o_f_null")
        colnames(tab_5)=c("statistic","series","Estimated_distribution")
        
        tab_lines=rbind(tab_3,tab_4,tab_5)
        
        maxim_entry_hists=max(tab_hists$series)
        maxim_entry_lines=max(tab_lines$series)
        
        
        ### Generate X axis tickmarks: Intercalate among x axis tickmarks one at p*
     
        if(10^(-2)<p_star)
        {
            p_star_write=format(round(p_star, 2))
            if(p_star==1) {p_star_write=1}
        }else{
            p_star_write=format(p_star,scientific=TRUE)
            if(p_star==0) {p_star_write=0}
        }
        
        if(10^(-2)<rho_o_scaling)
        {
            rho_o_write=round(rho_o_scaling,2)
            if(rho_o_scaling==1) {rho_o_write=1}

        }else{
            rho_o_write=format(rho_o_scaling,digits=2,scientific=TRUE)
            if(rho_o_scaling==0) {rho_o_write=0}
        }
        
        if(10^(-2)<pi_o)
        {
            pi_o_write=round(pi_o,2)
            if(pi_o==1) {pi_o_write=1}
        }else{
            pi_o_write=format(pi_o,digits=2,scientific=TRUE)
            if(pi_o==0) {pi_o_write=0}
        }
        
        
        tickmarks_seq=seq(min_stat,end,(end-min_stat)/5)
        
        index_p_star_ticks=which.min(abs(tickmarks_seq-p_star))
        tickmarks_length=length(tickmarks_seq)
        
        if(tickmarks_seq[index_p_star_ticks]>p_star)
        {
            if(min(abs(tickmarks_seq-p_star))>0.09*(end-min_stat))
            {
                tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star,tickmarks_seq[index_p_star_ticks:length(tickmarks_seq)])
                tickmarks_length=tickmarks_length+1
            }else{
                if(index_p_star_ticks<length(tickmarks_seq))
                {
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star,tickmarks_seq[index_p_star_ticks+1:length(tickmarks_seq)])
                }else{
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star)
                }
            }
        }else{
            if(min(abs(tickmarks_seq-p_star))>0.09*(end-min_stat))
            {
                if(index_p_star_ticks<length(tickmarks_seq))
                {
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks],p_star,tickmarks_seq[(index_p_star_ticks+1):length(tickmarks_seq)])
                    tickmarks_length=tickmarks_length+1
                    index_p_star_ticks=index_p_star_ticks+1
                }else{
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks],p_star)
                    tickmarks_length=tickmarks_length+1
                    index_p_star_ticks=index_p_star_ticks+1
                }
            }else{
                if(index_p_star_ticks<length(tickmarks_seq))
                {
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star,tickmarks_seq[index_p_star_ticks+1:length(tickmarks_seq)])
                }else{
                    tickmarks_seq=c(tickmarks_seq[1:(index_p_star_ticks-1)],p_star)
                }
            }
        }
        
        if(end>0.1)
        {
            labels_tickmarks=round(tickmarks_seq,2)
        }else{
            labels_tickmarks=format(tickmarks_seq,digits=3,scientific=TRUE)
        }
       
        labels_tickmarks[index_p_star_ticks]=paste0(p_star_write)
        #labels_tickmarks[index_p_star_ticks]=paste0("p*=\n",p_star_write)
        
        if(all(flag_grenander_choice==c(1,4))|all(flag_grenander_choice==c(4,1)))
        {
            x_annotation=(end-min_stat)*0.7
        }else if(all(flag_grenander_choice==c(3,4))|all(flag_grenander_choice==c(4,3))|all(flag_grenander_choice==c(4,4))){
            x_annotation=(end-min_stat)*0.3
        }else{
            x_annotation=(end-min_stat)*0.5}
        
        p1<-ggplot() +
        geom_area(data=tab_hists, aes(x=statistic,y=series,fill=Data),alpha=0.5,position="identity", na.rm=TRUE)+
        geom_line(data=tab_lines,aes(x=statistic,y=series,colour=Estimated_distribution),size=0.3,position="identity", na.rm=TRUE)+
        scale_fill_manual(values=c("firebrick2","dodgerblue1"),labels=c("Full"="Actual tests","Permuted"="Permut. tests"))+
        scale_colour_manual(values=c("firebrick","dodgerblue1","grey50"),labels=c("f_full"=expression(paste("Actual tests ",f)),"f_perm"=expression(paste("Permut. tests ",f[{o}])),"pi_o_f_null"=expression(paste("Estimated null ",pi[{o}],"·",f[{null}]))))+
        ylab("Local densities")+
        xlab("P-value")+
        guides(colour=guide_legend(title="Local densities",label.hjust = 0, label.vjust = 0, override.aes = list(size=1.5)))+
        guides(fill=FALSE)+
        #guides(fill=guide_legend(title="Local densities",override.aes = list(size=3)))+
        theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.x  = element_text(size=10),
        axis.title.y  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))+
        scale_x_continuous(breaks=tickmarks_seq,labels=labels_tickmarks,limits = c(min_stat*0.999,end))+
        geom_vline(xintercept=p_star, colour="grey", linetype = "longdash", na.rm=TRUE)+
        annotate("text", x = x_annotation, y=maxim*0.9, label =paste("f", "== paste(pi[{o}],f[{null}],+(1-pi[{o}]),f[{alt}])"),parse=TRUE)
        
        ####################
        
        if(full_data_is_file){
            x=data.frame(statistic=grid_unif,F=F_unif,F_o=F_o_unif)
            y=rho_fit_table
            colnames(y)=c("statistic","rho_hat","pi_fit")
            x=x[which(x$statistic<end),]
        }else{
            x=data.frame(statistic=grid,F=F,F_o=F_o)
            y=rho_fit_table
            colnames(y)=c("statistic","rho_hat","pi_fit")
            x=x[which(x$statistic<end),]
        }
        #ECDFs
        x=melt(x,id="statistic")
        #splines_data
        y=melt(y,id="statistic")
        
        #pi_o_fit_data:just the dots
        z=y[which(y$variable=="rho_hat"),]
        fact=100
        if(nrow(z)>fact)
        {
            set_step=as.integer(nrow(z)/fact)
            z=z[seq(1,nrow(z),set_step),]
        }
        
        #pi_o_fit_data:just the fit
        zz=y[which(y$variable=="pi_fit"),]
        
        #x=rbind(x,zz)
        x$variable=factor(x$variable)
        
        maxim_cumulative=max(x$value[which(x$statistic<=end)])
        
        tickmarks_seq_y_distributions=seq(0,maxim_cumulative*1.01,maxim_cumulative/5)
        labels_tickmarks_seq_y_distributions=round(tickmarks_seq_y_distributions,2)
        
        p2<-ggplot() +
        geom_line(data=x, aes(x=statistic,y=value,colour=variable),size=0.7, na.rm=TRUE)+
        scale_colour_manual(values=c("dodgerblue1","firebrick","mediumpurple4"),labels=c("F"=expression(F),"F_o"=expression(F[{o}]),"pi_fit"=expression(paste(tilde(pi)[{o}],"  fit"))))+
        ylab("Cumulative distributions")+
        xlab("P-value")+
        guides(colour=guide_legend(title="CDFs",label.hjust = 0, label.vjust = 0,override.aes = list(size=1.5),order = 1),
        fill=guide_legend(title="",override.aes = list(size=2),order = 2))+
        theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))+
        scale_x_continuous(breaks=tickmarks_seq,labels=labels_tickmarks,limits = c(min_stat*0.999,end))+
        scale_y_continuous(breaks=tickmarks_seq_y_distributions,labels=labels_tickmarks_seq_y_distributions,limits = c(-0.0001,maxim_cumulative*1.01))+
        geom_vline(xintercept=p_star, colour="grey", linetype = "longdash", na.rm=TRUE)
        
        ##########################################################################
        
        rm(tickmarks_seq)
        #end_pi_o_x_axis=(max(p_star,end,as.numeric(p_star_exhaustive_sweep==FALSE)*max_stat))
        end_pi_o_x_axis=max_stat
        tickmarks_seq=seq(min_stat,end_pi_o_x_axis,(end_pi_o_x_axis-min_stat)/5)
        index_p_star_ticks=which.min(abs(tickmarks_seq-p_star))
        tickmarks_length=length(tickmarks_seq)
        
        if(tickmarks_seq[index_p_star_ticks]>p_star)
        {
            if(min(abs(tickmarks_seq-p_star))>0.1*(end_pi_o_x_axis-min_stat))
            {
                tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star,tickmarks_seq[index_p_star_ticks:length(tickmarks_seq)])
                tickmarks_length=tickmarks_length+1
            }else{
                if(index_p_star_ticks<length(tickmarks_seq))
                {
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star,tickmarks_seq[index_p_star_ticks+1:length(tickmarks_seq)])
                }else{
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star)
                }
            }
        }else{
            if(min(abs(tickmarks_seq-p_star))>0.1*(end_pi_o_x_axis-min_stat))
            {
                if(index_p_star_ticks<length(tickmarks_seq))
                {
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks],p_star,tickmarks_seq[(index_p_star_ticks+1):length(tickmarks_seq)])
                    tickmarks_length=tickmarks_length+1
                    index_p_star_ticks=index_p_star_ticks+1
                }else{
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks],p_star)
                    tickmarks_length=tickmarks_length+1
                    index_p_star_ticks=index_p_star_ticks+1
                }
            }else{
                if(index_p_star_ticks<length(tickmarks_seq))
                {
                    tickmarks_seq=c(tickmarks_seq[1:index_p_star_ticks-1],p_star,tickmarks_seq[index_p_star_ticks+1:length(tickmarks_seq)])
                }else{
                    tickmarks_seq=c(tickmarks_seq[1:(index_p_star_ticks-1)],p_star)
                }
            }
        }
        
        if(end_pi_o_x_axis>0.1)
        {
            labels_tickmarks=round(tickmarks_seq,2)
        }else{
            labels_tickmarks=format(tickmarks_seq,digits=3,scientific=TRUE)
        }
        labels_tickmarks[index_p_star_ticks]=paste0(p_star_write)
        
        tickmarks_seq_y=c(0,0.25,0.50,0.75,1.00)
        index_pi_o_ticks=which.min(abs(tickmarks_seq_y-rho_o_scaling))
        
        if(tickmarks_seq_y[index_pi_o_ticks]>rho_o_scaling)
        {
            if(min(abs(tickmarks_seq_y-rho_o_scaling))>0.1)
            {
                tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],rho_o_scaling,tickmarks_seq_y[index_pi_o_ticks:length(tickmarks_seq_y)])
            }else{
                if(index_pi_o_ticks<length(tickmarks_seq_y)){
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],rho_o_scaling,tickmarks_seq_y[index_pi_o_ticks+1:length(tickmarks_seq_y)])
                }else{
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],rho_o_scaling)
                }
            }
        }else{
            if(min(abs(tickmarks_seq_y-rho_o_scaling))>0.1)
            {
                if(index_pi_o_ticks<length(tickmarks_seq_y)){
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks],rho_o_scaling,tickmarks_seq_y[index_pi_o_ticks+1:length(tickmarks_seq_y)])
                }else{
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks],rho_o_scaling)
                }
                
            }else{
                if(index_pi_o_ticks<length(tickmarks_seq_y)){
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],rho_o_scaling,tickmarks_seq_y[index_pi_o_ticks+1:length(tickmarks_seq_y)])
                }else{
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],rho_o_scaling)
                }
            }
        }
        labels_tickmarks_seq_y=round(tickmarks_seq_y,2)
        index_pi_o_ticks=which.min(abs(tickmarks_seq_y-rho_o_scaling))
        eval(parse(text=paste0("labels_tickmarks_seq_y[index_pi_o_ticks]=expression(rho[{o}]==",rho_o_write,")")))
        
        if(rho_o_scaling<0.5){
            x_legend=0.55
            y_legend=0.95
        }else if(p_star>0.5)
        {
            x_legend=0.45
            y_legend=0.4
        }else{
            x_legend=0.55
            y_legend=0.4
        }
        
        if(p_star_exhaustive_sweep==TRUE){
            
            pi_o_plot=data.frame(statistic=rho_o_scaling_loop_table$p_star_candidate,variable="pi_o_scalingloop",value=rho_o_scaling_loop_table$rho_o_scaling_candidate)

            z=rbind(z,pi_o_plot)
            lbl_o=paste0('paste("Scaling p-value ",p^symbol("\052")==arg.min(rho[{o}](p)),"=',p_star_write,'")')
            #Permuted distribution scaling factor ",rho[{o}]==min(rho[{o}](p)))'

            lbl=paste0('paste("Permuted distribution scaling factor ",rho[{o}]==min(rho[{o}](p)),"=',rho_o_write,'")')
            lbl_2=paste0('paste("Fraction of true null tests ",pi[{o}]==rho[{o}]*F[{o}](p^symbol("\052"))+(1-F(p^symbol("\052"))),"=',pi_o_write,'")')
            
            #lbl_o='paste("Scaling p-value ",p^symbol("\052")==arg.min(rho[{o}](p)))'
            #Permuted distribution scaling factor ",rho[{o}]==min(rho[{o}](p)))'
            
            #lbl='paste("Permuted distribution scaling factor ",rho[{o}]==min(rho[{o}](p)))'
            #lbl_2=paste0('paste("Fraction of true null tests ",pi[{o}]==rho[{o}]*F[{o}](p^symbol("\052"))+(1-F(p^symbol("\052"))),"=',pi_o,'")')
            
            
            p3<-ggplot() +
            geom_line(data=zz, aes(x=statistic,y=value),size=0.7,colour="dodgerblue1", na.rm=TRUE)+
            geom_point(data=z, aes(x=statistic,y=value,fill=variable),shape=21,size=0.8, na.rm=TRUE)+
            scale_colour_continuous(guide = FALSE)+
            scale_fill_manual(values=c("dodgerblue1","firebrick"),labels=c("rho_hat"=expression(paste(tilde(rho)[{o}](p,p^symbol("\052"))==frac(F(p^symbol("\052"))-F(p),F[{o}](p^symbol("\052"))-F[{o}](p)))),"rho_o_scaling_vector"=expression(paste(rho[{o}],"(p)=",lim(tilde(rho)[{o}](p*symbol("\242"),p), p*symbol("\242") %->% p)))))+
            ylab(expression(paste("Scaling factors  ",tilde(rho)[{o}]," , ",rho[{o}])))+
            xlab("P-value")+
            guides(fill=guide_legend(title="",override.aes = list(size=2),label.hjust = 0, order = 2))+
            theme(axis.text.y   = element_text(size=10),
            axis.text.x   = element_text(size=10),
            axis.title.y  = element_text(size=10),
            axis.title.x  = element_text(size=10),
            #legend.position = c(x_legend, y_legend),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.title=element_text(size=12),
            legend.key.size = unit(3, "lines"),
            legend.text=element_text(size=12))+
            scale_x_continuous(breaks=tickmarks_seq,labels=labels_tickmarks,limits = c(min_stat*0.999,end_pi_o_x_axis))+
            scale_y_continuous(breaks=tickmarks_seq_y,labels=labels_tickmarks_seq_y,limits = c(-0.0001,1.0001))+
            geom_vline(xintercept=p_star, colour="grey", linetype = "longdash", na.rm=TRUE)+
            geom_hline(yintercept=rho_o_scaling, colour="grey", linetype = "longdash", na.rm=TRUE)+
            annotate("text", x = x_legend, y=y_legend, label =lbl_o,size=3,parse=TRUE)+
            annotate("text", x = x_legend, y=y_legend-0.17, label =lbl,size=3,parse=TRUE)+
            annotate("text", x = x_legend, y=y_legend-0.3, label =lbl_2,size=3,parse=TRUE)
        }else{
            
            lbl_o=paste0('paste("Scaling p value ",p^symbol("\052")==arg.min(f(p)/f[{o}](p)),"=',p_star_write,'")')
            lbl=paste0('paste("Permuted distribution scaling factor ",rho[{o}](p^symbol("\052"))==lim(tilde(rho)[{o}](p,p^symbol("\052")), p %->% p^symbol("\052")),"=',rho_o_write,'")')
            lbl_2=paste0('paste("Fraction of true null tests ",pi[{o}](p^symbol("\052"))==rho[{o}](p^symbol("\052"))*F[{o}](p^symbol("\052"))+(1-F(p^symbol("\052"))),"=',pi_o_write,'")')
            
            p3<-ggplot() +
            geom_line(data=zz, aes(x=statistic,y=value),size=0.7,colour="dodgerblue1", na.rm=TRUE)+
            geom_point(data=z, aes(x=statistic,y=value,fill=variable),shape=21,color="grey95",size=0.8, na.rm=TRUE)+
            scale_colour_continuous(guide = FALSE)+
            scale_fill_manual(values=c("dodgerblue1"),labels=c("rho_hat"=expression(paste(tilde(rho)[{o}](p,p^symbol("\052"))==frac(F(p^symbol("\052"))-F(p),F[{o}](p^symbol("\052"))-F[{o}](p))))))+
            ylab(expression(paste("Scaling factors  ",tilde(rho)[{o}]," , ",rho[{o}])))+
            xlab("P-value")+
            guides(fill=guide_legend(title="",override.aes = list(size=2),label.hjust = 0, order = 2))+
            theme(axis.text.y   = element_text(size=10),
            axis.text.x   = element_text(size=10),
            axis.title.y  = element_text(size=10),
            axis.title.x  = element_text(size=10),
            #legend.position = c(x_legend, y_legend),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.title=element_text(size=12),
            legend.key.size = unit(3, "lines"),
            legend.text=element_text(size=12))+
            scale_x_continuous(breaks=tickmarks_seq,labels=labels_tickmarks,limits = c(min_stat*0.999,end_pi_o_x_axis))+
            scale_y_continuous(breaks=tickmarks_seq_y,labels=labels_tickmarks_seq_y,limits = c(-0.0001,1.0001))+
            geom_vline(xintercept=p_star, colour="grey", linetype = "longdash", na.rm=TRUE)+
            geom_hline(yintercept=rho_o_scaling, colour="grey", linetype = "longdash", na.rm=TRUE)+
            annotate("text", x = x_legend, y=y_legend, label =lbl_o,size=3,parse=TRUE)+
            annotate("text", x = x_legend, y=y_legend-0.17, label =lbl,size=3,parse=TRUE)+
            annotate("text", x = x_legend, y=y_legend-0.3, label =lbl_2,size=3,parse=TRUE)
            
            
            
        }
        ############################
        ############################
        
        cut=which(fdrs$statistic>=end*0.999999)[1]
        fdrs=fdrs[1:cut,c(1,(ncol(fdrs):2))]
        tab=tab[1:cut,]
        
        labels_fdr=labels_fdr[length(labels_fdr):1]
        
        colors_number=ncol(fdrs)-2
        palette_set_1=c("dodgerblue1","firebrick","mediumpurple4","#4DAF4A","#FF7F00","#FFFF33","#A65628","#F781BF")
        palette=palette_set_1[colors_number:1]
        
        fdrs=fdrs[which(fdrs$statistic<end),]
        fdrs=fdrs[,c(2:ncol(fdrs))]
        fdrs=melt(fdrs,id="quantile")
        
        maxim_fdrs=max(fdrs$value[which(fdrs$variable %in% c("Fdr_MV","Fdr_BH_perm","Fdr_BH_unif","Fdr_ST","Fdr"))])
        
        
        threshold_tick=significance_threshold
        tickmarks_seq_y=c(0,0.25,0.50,0.75,1.00)
        index_pi_o_ticks=which.min(abs(tickmarks_seq_y-threshold_tick))
        
        if(tickmarks_seq_y[index_pi_o_ticks]>threshold_tick)
        {
            if(min(abs(tickmarks_seq_y-threshold_tick))>0.1)
            {
                tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],threshold_tick,tickmarks_seq_y[index_pi_o_ticks:length(tickmarks_seq_y)])
            }else{
                if(index_pi_o_ticks<length(tickmarks_seq_y)){
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],threshold_tick,tickmarks_seq_y[index_pi_o_ticks+1:length(tickmarks_seq_y)])
                }else{
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],threshold_tick)
                }
            }
        }else{
            if(min(abs(tickmarks_seq_y-threshold_tick))>0.1)
            {
                if(index_pi_o_ticks<length(tickmarks_seq_y)){
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks],threshold_tick,tickmarks_seq_y[index_pi_o_ticks+1:length(tickmarks_seq_y)])
                }else{
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks],threshold_tick)
                }
                
            }else{
                if(index_pi_o_ticks<length(tickmarks_seq_y)){
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],threshold_tick,tickmarks_seq_y[index_pi_o_ticks+1:length(tickmarks_seq_y)])
                }else{
                    tickmarks_seq_y=c(tickmarks_seq_y[1:index_pi_o_ticks-1],threshold_tick)
                }
            }
        }
        
        if(maxim_fdrs>0.1)
        {
            labels_tickmarks_seq_y=round(tickmarks_seq_y,2)
        }else{
            labels_tickmarks_seq_y=format(tickmarks_seq_y,digits=3,scientific=TRUE)
        }
        
        
        
        #tickmarks_seq_y_fdrs=seq(0,maxim_fdrs*1.01,maxim_fdrs/5)
        #labels_tickmarks_seq_y_fdrs=round(tickmarks_seq_y_fdrs,2)
        p4<-ggplot() +
        geom_line(data=fdrs, aes(x=quantile,y=value,colour=variable),size=0.7, na.rm=TRUE)+
        scale_colour_manual(values=palette,labels=labels_fdr)+
        ylab("Fdr & Fndr")+
        guides(colour=guide_legend(title="Methods",label.hjust = 0, label.vjust = 0,override.aes = list(size=1.5),order = 1))+
        ylab("Fdr & Fndr")+
        xlab("Number of tests")+
        theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))+
        scale_y_continuous(breaks=labels_tickmarks_seq_y,labels=labels_tickmarks_seq_y,limits = c(-0.0001,maxim_fdrs*1.01))+
        geom_vline(xintercept=p_star, colour="grey", linetype = "longdash", na.rm=TRUE)+
        geom_hline(yintercept=significance_threshold, colour="grey", linetype = "longdash", na.rm=TRUE)
        #annotate("text", x = (max(fdrs$quantile))*0.9, y=significance_threshold+0.05*maxim_fdrs, label =paste0(as.integer(significance_threshold*100,0),"% FDR"))
        
        pdf(paste0(output_name,"/fdrs.pdf"),height=9,width=8,onefile=FALSE)
        print(perm_fdr_plot(list(p1,p2,p3,p4)))
        graphics.off()
    }
    
    if(print_messages)
    {
        message("#############################################################################################################")
        message("###################                                  Done!                                 ##################")
        message("#############################################################################################################")
    }
    
    return(output)
}
