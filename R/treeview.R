
#' Generate interactive tree visualisations and scatter plots for illustrating scannint statistics. 
#' 
#' This will produce a set of html widgets which will highlight by colour and tooltips statistics such as growth rate and molecular clock outliers. 
#' 
#' @param e0 Path to the scanner environment produced by \code{tfpscan}. Alternatively can pass the environment directly. 
#' @param branch_cols A character vector of statistics for which outputs should be produced. The logistic growth rate plot will always be produced. 
#' @param mutations A character vector of mutations which will be illustrated in a heatmap
#' @param lineages A set of lineage names which will be used to subdivide outputs in scatter plots. 
#' @param output_dir Outputs will be saved in this directory. Will create the directory if it does not exist. 
#' @return A ggtree plot
#' @export
treeview <- function( e0
, branch_cols = c('logistic_growth_rate', 'clock_outlier')
, mutations = c( 'S:A222V', 'S:Y145H', 'N:Q9L', 'S:E484K')
, lineages = c( 'AY\\.9' , 'AY\\.43', 'AY\\.4\\.2')
, output_dir = 'treeview' 
)
{
	library( ggtree ) 
	library( ape ) 
	library( lubridate ) 
	library( glue ) 
	library( ggplot2 )
	library( ggiraph )
	
	dir.create( output_dir, showWarnings = FALSE )
	
	# load env 
	if ( is.character( e0 ) ){
		e0 = readRDS(e0) 
	}
	sc0 <- e0$Y 
	cmuts = lapply( 1:nrow(sc0), function(i){
		list( 
		defining = strsplit( sc0$all_mutations[i], split = '\\|' )[[1]]
		, all = strsplit( sc0$defining_mutations[i], split = '\\|' )[[1]] 
		)
	})
	names( cmuts )<- sc0$cluster_id
	tr1 <- e0$tre
	
	stopifnot( all( branch_cols %in% colnames( e0$Y )) )
	
	# pick one tip from each cluster 
	sc0$representative = NA 
	for ( i in 1:nrow( sc0 )){
		if ( sc0$external_cluster[i])
		{
			tu = strsplit( sc0$tip[i], split = '\\|' )[[1]]
			.sts <- e0$sts[ tu ] 
			sc0$representative[i] <-  tu[ which.max( .sts  )  ] # most recent sample in cluster 
		}
	}
	tr2 = keep.tip( tr1, na.omit( sc0$representative ) ) 
	
	# subtrees 
	message( 'Computing sub-trees' )
	stres2 <- subtrees( tr2, wait = TRUE ) 
	Nstres2 <- sapply( stres2, Ntip )
	
	# for each cluster find the node in tr2 representing it
	## internal nodes first 
	for ( i in seq_along( stres2 ) ){
		inode = i + Ntip( tr2 )
		uv = na.omit( sc0$node_number[ match( stres2[[i]]$tip.label , sc0$representative ) ]  ) 
		shared_anc = Reduce( intersect, e0$ancestors[uv] )
		shared_anc2 = setdiff( intersect ( shared_anc , sc0$node_number ), uv )
		if ( length( shared_anc2 ) > 0 ){
			a = shared_anc2[ which.min( e0$ndesc[shared_anc2] ) ]
			sc0$tr2mrca[ sc0$node_number == a ] <- inode 
		}
	}
	## tips (takes precedence if overlap in tr2mrca ) 
	i <- which( sc0$representative %in% tr2$tip.label )
	sc0$tr2mrca[i] <- match( sc0$representative[i], tr2$tip.label )
	
	# tree data frames 
	## tips 
	sc2 <- sc0[ !is.na( sc0$tr2mrca ) , ]
	sc2$date_range <- sapply( 1:nrow( sc2 ) , function(i) glue( '{sc2$least_recent_tip[i]} -> {sc2$most_recent_tip[i]}') )
	tdvars = unique( c(branch_cols, 'logistic_growth_rate', 'clock_outlier', 'cluster_size', 'date_range', 'cluster_id', 'region_summary', 'cocirc_lineage_summary', 'lineage' , 'tr2mrca')  )
	td0 <- sc2[ sc2$tr2mrca <= Ntip( tr2 ) 
	 , tdvars  
	]
	td0$lineages = td0$lineage 
	td0$cocirc_summary = td0$cocirc_lineage_summary 
	td0$node = td0$tr2mrca 
	td0$internal = 'N' 
	##internal
	td1 <- sc2[ sc2$tr2mrca > Ntip( tr2 ) 
	, tdvars
	]
	if ( nrow( td1 ) > 0 ){
		td1$lineages = td1$lineage 
		td1$cocirc_summary = td1$cocirc_lineage_summary 
		td1$node = td1$tr2mrca 
		td1$internal = 'Y' 
		td1$cluster_size <- 0 
		x = setdiff( (Ntip(tr2)+1):(Ntip(tr2) + Nnode(tr2)) ,   td1$node ) # make sure every node represented 
		td1 = merge( td1, data.frame( node = x ), all =TRUE )
		td <- rbind( td0, td1 )
	} else{
		td = td0 
	}
	td = td[ order( td$node ), ] # important 

	# rescale clock ?
	td$clock_outlier <- scale(  td$clock_outlier )/2
	
	# interpolate missing values  &  repair cluster sizes 
	td$logistic_growth_rate[ (td$node <= Ntip( tr2 )) & (is.na(td$logistic_growth_rate)) ] <- 0 
	td$clock_outlier[ (td$node <= Ntip( tr2 )) & (is.na(td$clock_outlier)) ] <- 0 
	for (ie in postorder( tr2 )){
		a = tr2$edge[ie,1]
		u = tr2$edge[ie,2]
		td$cluster_size[a] <- td$cluster_size[a] + td$cluster_size[u]
		for ( vn in branch_cols ) {
			if ( is.na( td[[ vn ]][ a] ) ){
				td[[ vn ]] [ a ]  <- td[[vn ]] [ u ] 
			}
		}
	}
	
	
	# cols for continuous stats  
	cols = rev( c("red", 'orange', 'green', 'cyan', 'blue') )
	
	# lineages for clade labels 
	td$lineages1 <- sapply( strsplit( td$lineages, split = '\\|') , '[', 1 )
	sc0$lineage1 <- sapply( strsplit( sc0$lineage, split = '\\|') , '[', 1 )
	
	tablin <- table(  td$lineages1[ !(sc0$lineage1 %in% c('None', 'B.1')) ]  )
	tablin <- tablin [order( tablin ) ]
	lins <- names( tablin )
	
	# find a good internal node to represent the mrca of each lineage 
	lin_nodes = c() 
	lin_node_names <- c() 
	for (lin in lins){
		whichrep = na.omit( sc0$representative[ sc0$lineage1==lin ]  ) 
		if ( length( whichrep) == 1){
			res = which( tr2$tip.label == whichrep )
		}else{
			res = getMRCA( tr2, whichrep )
		}
		if ( !is.null(res)){
			lin_nodes <- c( lin_nodes, res )
			lin_node_names <- c( lin_node_names, lin )
		}
		#lin_node_names <- lin_node_names[!duplicated(lin_nodes) ]
		#lin_nodes <- lin_nodes[!duplicated(lin_nodes) ]
		
		res
	}
	
	# make and save the tree view 
	.plot_tree <- function (vn , mut_regex = NULL , colour_limits = NULL )
	{
		if ( is.null( colour_limits )){
			colour_limits = range(td[[vn]] ) 
		}
		gtr1 = ggtree( dplyr::full_join( tr2, td , by = 'node') , aes_string( colour = vn) , ladderize = TRUE, right = TRUE, continuous = TRUE) 
		
		shapes = c( Y = '\U2B24', N = '\U25C4' )
		gtr1.1 <- gtr1 +
		 scale_color_gradientn( name = gsub(vn, patt = '_', rep = ' '), colours  = cols, limits = colour_limits, oob = scales::squish ) + 
		 geom_point( aes_string(color = vn, size = 'cluster_size', shape = 'as.factor(internal)'),  data = gtr1$data) + 
		 scale_shape_manual( name = NULL, labels = NULL, values = shapes ) +
		 scale_size(name = 'Cluster size', range = c(2, 16)) + 
		 ggtitle( glue( '{Sys.Date()}, colour: {vn}')  )  + theme(legend.position='top' )
		
		for( i in  seq_along(lins) ){
			if ( !is.na( lin_nodes[i] ) ) {
				gtr1.1 <- gtr1.1 + geom_cladelabel(node=lin_nodes[i], label = lin_node_names[i], offset = .00001, colour = 'black') 			
			}
		}
		
		ggsave( gtr1.1, file = glue('{output_dir}/tree-{vn}-{Sys.Date()}.svg') 
			,  height = max( 14, floor( Ntip( tr2 ) / 10 ) )
			, width = 16 
			, limitsize = FALSE )
		
		# make mouseover info  
		## standard meta data
		ttdfs = apply( gtr1.1$data, 1, FUN = function(x){
			z= as.list( x )
			lgr = as.numeric( z$logistic_growth_rate )
			y = with (z ,
				data.frame(
				  `Cluster ID` = glue( '#{cluster_id}') 
				  , `Cluster size` = cluster_size 
				  , `Date range` = date_range
				  , `Example sequence` = label 
				  , `Logistic growth` =   paste0( ifelse(lgr>0, '+', ''), round( lgr*100 ) , '%')
				  , `Mol clock outlier` = clock_outlier 
				  #, `Structure Z` = treestructure_z
				  , `Lineages` = lineages
				)
			)
			y = t(y) 
			colnames(y) <- '' 
			tryCatch( paste( knitr::kable( y , 'simple' ) , collapse = '\n' )
				, error = function(e)  paste( knitr::kable( y , 'markdown' ) , collapse = '\n' ) )
		})


		## table with geo composition 
		ttregtabs = gtr1.1$data$region_summary # 
		## cocirc
		ttcocirc = gtr1.1$data$cocirc_summary #
		## defining muts 
		.sort.muts <- function( muts ){
			if ( length( muts ) == 0 )
				return('') 
			pre <- sapply( strsplit( muts, split = ':') , '[', 1 )
			upres <- sort( unique( pre ))
			do.call( c, lapply( upres, function(.pre){
				.muts <- muts[ pre == .pre ]
				.muts1 <- sapply( strsplit( .muts, split = ':'), '[', 2 )
				sites <- regmatches( .muts1, regexpr(.muts1, patt = '[0-9]+' ) )
				o <- order( as.numeric( sites ))
				.muts[o]
			}))
		}
		ttdefmuts <- sapply( match( gtr1.1$data$cluster_id, sc0$cluster_id ) , function(isc0){
			if ( is.na( isc0 ))
				return( '' ) 
			paste( sep = '\n' , 'Cluster branch mutations:',
				gsub( 
					x = tryCatch(
						stringr::str_wrap( 
							paste( collapse = ' ' , .sort.muts( cmuts[[  as.character( sc0$node_number[isc0] ) ]]$defining )  )
							, width = 60
						)
						, error = function(e)browser() )
					, patt = ' '
					, rep = ', '
				)
				,'\n'
			)
		})
		ttallmuts <- sapply( match( gtr1.1$data$cluster_id, sc0$cluster_id ) , function(isc0){
			if ( is.na( isc0 ))
				return( '' ) 
			paste( sep = '\n' , 'All mutations:',
				gsub( 
					x = stringr::str_wrap( 
						paste( collapse = ' ' , .sort.muts(  cmuts[[  as.character( sc0$node_number[isc0] ) ]]$all ) )
						, width = 60
					)
					, patt = ' '
					, rep = ', '
				)
				,'\n'
			)
		})
		
		gtr1.1$data$defmuts = ttdefmuts 
		gtr1.1$data$allmuts = ttallmuts 
		if ( !is.null( mut_regex  )){
			#gtr1.1$data$label <- '' 
			for ( mre in mut_regex){
				i <- which( grepl(gtr1.1$data$allmuts, patt = mre ) ) 
				#gtr1.1$data$label[i] <- '*'
				gtr1.1$data[[ mre ]] <-  grepl(gtr1.1$data$allmuts, patt = mre )
			}
			#gtr1.1 <- gtr1.1 + geom_tiplab( size = 16, colour = 'red' )
		}
		
		genotype <- as.data.frame( gtr1.1$data[ gtr1.1$data$node <= Ntip(tr2) , c( 'label', mut_regex ) ] )
		rownames( genotype ) <- genotype$label 
		genotype <- genotype[ ,-1, drop = FALSE] 
		
		# make html widget 
		gtr1.1$data$mouseover = sapply( 1:length( ttdfs ), function(i){
			 paste0(  'Statistics:\n',  ttdfs[i],   '\n\nGeography:\n', ttregtabs[i], '\n\nCo-circulating with:\n', ttcocirc[i], '\n\n', ttdefmuts[i],'\n', ttallmuts[i],'\n' , collapse = '\n')
		})
		gtr1.1$data$colour_var <-gtr1.1$data[[ vn ]]
		gtr1.2 <- gheatmap( gtr1.1, genotype , width = .075 , offset=0.0005, colnames_angle=-90, colnames_position='top', colnames_offset_y=-6, legend_title = 'Genotype')
		gtr1.3 <- gtr1.2 +   
		ggiraph::geom_point_interactive(aes(x = x, y = y
			  , color = colour_var
			  , tooltip = mouseover
			  , data_id = node
			  , size = cluster_size+1
			  , shape = as.factor(internal)
			)
		) + scale_shape_manual( name = NULL, labels = NULL, values = shapes ) +  
		scale_size(name = 'Cluster size', range = c(2, 16)) + 
		scale_color_gradientn(  name = stringr::str_to_title( gsub(vn, patt = '_', rep = ' '))
			, colours  = cols
			, limits = colour_limits
			, oob = scales::squish ) + 
		theme( legend.position = 'top' )
		
		
		#~ font-family: "Lucida Console", "Courier New", monospace;
		tooltip_css <- "background-color:black;color:grey;padding:14px;border-radius:8px;font-family:\"Courier New\",monospace;"
		pgtr1.3  <- girafe(
		  ggobj = gtr1.3, width_svg = 15
		  , height_svg = max( 14,  floor( Ntip( tr2 ) / 10 ) )
		  ,  options = list(
						  opts_sizing(width = 0.8), 
						  opts_tooltip(css = tooltip_css, use_fill=FALSE)
						)
		)
		
		htmlwidgets::saveWidget( pgtr1.3 
		 , file = as.character( glue( '{output_dir}/tree-{vn}.html' ) )
		 , title = glue('SARS CoV 2 scan {Sys.Date()}'))
		file.copy( as.character( glue( '{output_dir}/tree-{vn}.html' ) ), as.character( glue( '{output_dir}/tree-{vn}-{Sys.Date()}.html' ) ) , overwrite = TRUE)
		
		gtr1.1
	}
	
	
	#' @param pldf the data element of a tree plot 
	.pl_cluster_sina <- function(pldf, varx = 'logistic_growth_rate'
	, mut_regexp  = 'S:A222V' , lineage_regexp = NULL 
	)
	{
		#~ 	pldf = pl$data 
		sc1 <- pldf [ pldf$isTip , ]
		
		sc1$varx = sc1[[ varx ]] 
		sc1$colour_var <- ''
		
		if ( !is.null( mut_regexp )){
			y = do.call( cbind, lapply( mut_regexp , function( xre ){
				ifelse( grepl( sc1$allmuts, patt = xre ),  xre, '' )
			}))
			ymut = apply( y, 1, function(x) paste( x, collapse = '.' ) )
		} else{
			ymut <- rep( '', nrow( sc1 ))
		}
		
		if ( !is.null( lineage_regexp )){
			y = do.call( cbind, lapply( lineage_regexp , function( xre ){
				ifelse( grepl( sc1$lineage, patt = xre ),  xre, '' )
			}))
			ylin = apply( y, 1, function(x) paste( x, collapse = '.' ) )
		} else {
			ylin = rep( '', nrow ( sc1 ))
		}
		
		sc1$mutation_lineage <- paste( ymut, ylin, sep = '_' )
		
		fx = as.factor( sc1$mutation_lineage )
		sc1$x <- as.numeric( fx )
		sc1$x <- rnorm( nrow( sc1 ) , sc1$x, .15 )
		sc1$y = sc1$logistic_growth_rate 
		
		p1 <- ggplot( sc1 ) + ggiraph::geom_point_interactive(aes( x = x, y = y, colour = mutation_lineage, size = cluster_size+1)
		  , tooltip = sc1$mouseover
		  , data_id = sc1$node
		  , alpha = .5
		  , data = sc1 
		) + xlab( 'Lineage and/or mutation' ) + ylab( 'Logistic growth rate' ) + geom_hline( aes( yintercept = 0 ) )
		
		tooltip_css <- "background-color:black;color:grey;padding:14px;border-radius:8px;font-family:\"Courier New\",monospace;"
		g1 = girafe(
		  ggobj = p1, width_svg =8
		  , height_svg = 8
		  , options = list(
						  opts_sizing(width = 0.8), 
						  opts_tooltip(css = tooltip_css, use_fill=FALSE)
						)
		)

		htmlwidgets::saveWidget( g1
		 , file = as.character( glue( '{output_dir}/sina-{varx}.html' ) )
		)
	}
	
	message( 'Generating figures' )
	pl = suppressWarnings( .plot_tree( 'logistic_growth_rate' , mut_regex = mutations , colour_limits = c(-.5,.5)) )
	pldf = pl$data 
	for ( vn in setdiff(  branch_cols, c('logistic_growth_rate') ) ){  
		suppressWarnings( .plot_tree( vn , mut_regex = mutations )  )
	}
	
	suppressWarnings( 
		.pl_cluster_sina(pldf 
		, mut_regexp  = mutations
		, lineage_regexp = lineages 
		)
	)
	
	invisible( pl )
}

#~ treeview( e0 = 'tfpscan-2021-11-25/scanner-env-2021-11-25.rds'
#~ , branch_cols = c('logistic_growth_rate')
#~ , mutations = c( 'S:Y145H', 'N:Q9L')
#~ , lineages = c( 'AY\\.9' , 'AY\\.43' )
#~ , output_dir = 'treeview' 
#~ )

