subroutine greedy(as_all, log_p_mat_ma, k_new, mRNA_all, n_genes, max_iter, nvars, optimize, p, &
                        sum_prob_tot, population, loops, expected_clusters, display, &
                        my_results_population, my_results_best, my_results_overallsum, &
						my_results_runtime_eachloop, my_results_param_idx, my_results_parameter, &
						clusterprobabilities,my_distance,cell_prob,opt_idx_a)


	implicit none
	
	
!----------------------------------------------------------------------------------------------!
!---------------------------------Declaring the variables--------------------------------------!
!----------------------------------------------------------------------------------------------!

	integer :: n_genes, max_iter, nvars, p, loops, expected_clusters, display
	!integer, dimension(:,:), allocatable :: as_all, mRNA_all, population, as_all_t, mRNA_all_t
	!
	integer, dimension(loops,nvars) :: as_all
	integer, dimension(nvars,n_genes) :: mRNA_all
	integer, dimension(loops, nvars) :: population
	real*8, dimension(1,loops) :: sum_prob_tot
	
	!
	!real*8, dimension(:,:), allocatable :: sum_prob_tot
	real*8 :: k_new(945,3), log_p_mat_ma(201,945), k_new_t(3,945), log_p_mat_ma_t(945,201)
	character(len=13) :: algorithm = 'greedy_cabsel'
	integer :: i, optimize
	
	integer my_results_population(loops,nvars)
	integer my_results_best(nvars)
	real*8    my_results_overallsum(loops)
	real*8    my_results_runtime_eachloop(loops)
	integer my_results_param_idx(expected_clusters,n_genes)
	real*8    my_results_parameter(expected_clusters,n_genes,3)
	
	real*8, dimension(nvars,expected_clusters) :: clusterprobabilities, my_distance, cell_prob
	
	integer opt_idx_a(expected_clusters,n_genes,loops)
		
	!allocate(as_all(loops,nvars))
	!allocate(mRNA_all(nvars,n_genes))
	!allocate(as_all_t(nvars,loops))
	!allocate(mRNA_all_t(n_genes,nvars))
	!allocate(population(loops,nvars))
	!allocate(sum_prob_tot(1,loops))
	
	population(:,:) = 0
	sum_prob_tot(1,:) = 0
	!print*, optimize
	!print*, as_all
	!print*, as_all
	
	!print*, mRNA_all
	!print*, log_p_mat_ma
	!print*, k_new
	
!----------------------------------------------------------------------------------------------!
!--------------------------------Reading from .txt files---------------------------------------!
!----------------------------------------------------------------------------------------------!
	
	!open(unit=10, file='/home/R/CALISTA-R-Fortran/as_all.txt')
    
    !read(10,*) as_all_t
    
    !as_all = transpose(as_all_t)
    
    !close(10)
    
    
    !open(unit=20, file='/home/R/CALISTA-R-Fortran/k_new.txt')
    
    !read(20,*) k_new_t
    
    !k_new = transpose(k_new_t)
    
    !close(20)
    
    
    !open(unit=30, file='/home/R/CALISTA-R-Fortran/log_p_mat_ma.txt')
    
    !read(30,*) log_p_mat_ma_t
    
    !log_p_mat_ma = transpose(log_p_mat_ma_t)
    
    !close(30)
    
    
    !open(unit=40, file='/home/R/CALISTA-R-Fortran/mRNA_all.txt')
    
    !read(40,*) mRNA_all_t
    
    !mRNA_all = transpose(mRNA_all_t)
    
    !close(40)
    
     
    call greedy_cabsel_f(as_all,log_p_mat_ma,k_new, &
    mRNA_all,n_genes,max_iter,nvars,optimize, &
    p,sum_prob_tot,population,loops, &
    expected_clusters,algorithm,display, &
    my_results_population, my_results_best, &
    my_results_overallsum, my_results_runtime_eachloop, &
    my_results_param_idx,my_results_parameter,clusterprobabilities,&
    my_distance,cell_prob)
    
    end subroutine greedy
!----------------------------------------------------------------------------------------------!
!-------------------------------Greedy Algorithm-----------------------------------------------!
!----------------------------------------------------------------------------------------------!
    
subroutine	greedy_cabsel_f(as_all,log_p_mat_ma,k_new, &
    mRNA_all,n_genes,max_iter,nvars,optimize, &
    p,sum_prob_tot,population,loops, &
    expected_clusters,algorithm,display, &
    my_results_population, my_results_best, &
    my_results_overallsum, my_results_runtime_eachloop, &
    my_results_param_idx,my_results_parameter,clusterprobabilities,&
    my_distance,cell_prob)
		
		integer :: n_genes, max_iter, nvars, p, loops, expected_clusters, display
		integer :: as_all(loops,nvars), mRNA_all(nvars,n_genes), population(loops,nvars)
		real*8 :: k_new(945,3), log_p_mat_ma(201,945), X(945,201)
		real*8 :: sum_prob_tot(1,loops)
		character(len=13) :: algorithm
		integer :: my_break, break_iter, check
		integer :: optimize, indx_best,best_as(nvars)
		real*8 startiteration(loops),enditeration(loops)
		
		integer q,jj,iterations,info,as(nvars), max_mRNA_counts, num_cells_mRNA(201,n_genes), comp(nvars)
		
		! as_sorted contains both sorted_as (col1) and sorted_idx (col2)
		integer sorted_as(nvars),sorted_idx(nvars),as_sorted(nvars,2),data_sorted(nvars,n_genes)
		
		integer clusters(expected_clusters), last_idx(expected_clusters), num_clusters, bounds(expected_clusters,2),idx,maxrel
		
		integer opt_idx_clusters(expected_clusters,n_genes), clust, idx_max_L(n_genes), arrayscal(1), idx_max_cell_prob(nvars)
		
		integer, dimension(:,:), allocatable :: cells_in_each_cluster
		
		integer, dimension(:), allocatable :: aaa
		real*8, dimension(:,:), allocatable :: bbb
		
		real*8 Y(201,n_genes), Z(945,n_genes), max_L(n_genes), Z_temp(945), logP_each_gene_each_clust(expected_clusters,n_genes)
		
		real*8 cell_prob(nvars,expected_clusters), my_distance(nvars,expected_clusters)
		real*8 cell_prob_new(expected_clusters)
		real*8 opt_param_each_gene(201), sum_prob, sum_prob_max(nvars)
		
		real*8 clusterprobabilities(nvars,expected_clusters), ccc(nvars,expected_clusters)
		real*8 log_clusterprobabilities(nvars,expected_clusters)
		
		real*8    my_results_clusterprobabilities(nvars,expected_clusters,loops)
		real*8    my_results_distance(nvars,expected_clusters,loops)
		real*8    my_results_cell_prob(nvars,expected_clusters,loops)
		real*8    my_results_sum_prob_tot(1,loops)
		real*8    my_results_overallsum(loops)
		real*8    my_results_parameter(expected_clusters,n_genes,3)
		real*8    my_results_runtime_eachloop(loops)
		integer my_results_population(loops,nvars)
		integer my_results_best(nvars)
		integer my_results_param_idx(expected_clusters,n_genes)
		integer opt_idx_a(expected_clusters,n_genes,loops)
		integer opt_id(expected_clusters,n_genes)
		real*8    optpar(expected_clusters,n_genes,3)
		
		!print*, "Welcome2 to Fortran!"
		
		!----------------------------------------------------------------------------------------------!
!-----------------------------First Greedy run (optimize = 1)----------------------------------!
!----------------------------------------------------------------------------------------------!
		
		if (optimize == 1) then
		
		do jj = 1,loops
		call cpu_time(startiteration(jj))
		!jj = 1
			as(:) = as_all(jj,:)
			max_mRNA_counts = 200
			num_cells_mRNA(:,:) = 0
			!my_break = 1
			!check = 1
				do iterations = 1,max_iter
				my_break = 1
			    check = 1
		
					call mysort(as,nvars,sorted_as,sorted_idx,as_sorted)
				
					do i = 1,nvars
					data_sorted(i,:) = mRNA_all(sorted_idx(i),:)
					end do
					
					do i = 1,expected_clusters
					clusters(i) = i
					end do
					
					j = 0
					
					do i = 1,nvars
						if (sorted_as(i)<sorted_as(i+1)) then
						j = j+1
						last_idx(j) = i
						end if
						last_idx(expected_clusters) = nvars
					end do
					
					num_clusters = size(clusters)
					
					bounds(:,2) = last_idx
					bounds(1,1) = 1
					
					do i = 2,num_clusters
						bounds(i,1) = bounds(i-1,2) + 1
					end do
					
					X = transpose(log_p_mat_ma)
					
					opt_idx_clusters(:,:) = 0
					logP_each_gene_each_clust(:,:) = 0
					
					do clust = 1,num_clusters
				
					IF( ALLOCATED(cells_in_each_cluster) )  DEALLOCATE( cells_in_each_cluster ) 						
					allocate(cells_in_each_cluster((bounds(clust,2)-bounds(clust,1)+1),n_genes))
				
					cells_in_each_cluster = data_sorted(bounds(clust,1):bounds(clust,2),:)	
					
					if (bounds(clust,1) == bounds(clust,2)) then
						print*, "Error: single cell in one cluster..."
					else

						call myhist(expected_clusters,n_genes,num_cells_mRNA,max_mRNA_counts,bounds,clust,cells_in_each_cluster)
					
					end if		
					

					 Y = num_cells_mRNA

					 Y(:,:) = Y(:,:) / sum(Y(:,1))
					 
					 !call multmat(X,Y,Z, n_genes, max_mRNA_counts)

					 Z = matmul(X,Y)

					 do i = 1,n_genes
						Z_temp = Z(:,i)
						max_L(i) = maxval(Z_temp)
						!arrayscal is defined to save the array output, we then convert rank 1 array to scalar...
						arrayscal = maxloc(Z_temp)
						idx_max_L(i) = arrayscal(1)
					 end do

					 opt_idx_clusters(clust,:) = idx_max_L
					 logP_each_gene_each_clust(clust,:) = max_L
					 
					 end do
					
					 idx_max_cell_prob(:) = 0
					 cell_prob(:,:) = 0.0
					 
					 my_distance = cell_prob	
					 					 
					 do i = 1,nvars
					 cell_prob_new(:) = 0.0
						do clust = 1,num_clusters
							do j = 1,n_genes
								opt_param_each_gene = log_p_mat_ma(:,opt_idx_clusters(clust,j))
								cell_prob_new(clust) = cell_prob_new(clust) + opt_param_each_gene(mRNA_all(i,j)+1)
							end do
						end do
						
						sum_prob = maxval(cell_prob_new)
						arrayscal = maxloc(cell_prob_new)
						idx_max_cell_prob(i) = arrayscal(1)
						sum_prob_tot(1,jj) = sum_prob_tot(1,jj) + sum_prob
						!maxrel = 1
						
						if (check == 1) then
							if (as(i) /= idx_max_cell_prob(i)) then
								my_break = 0
								check =0
							end if
						end if	
										
					 end do
					 
					as = idx_max_cell_prob
					population(jj,:) = as			 
					if (my_break == 1) exit
				
				!end of iter	
				end do	 
					
				opt_idx_a(:,:,jj) = opt_idx_clusters
				
				call cpu_time(enditeration(jj))
				
				!print*, jj
				
				!end of main loop (50 loops)
				end do
				
				!if (optimize==1) then
				
					if (maxval(population(1,:)) /= maxval(population(2,:))) then
						my_results_sum_prob_tot = sum_prob_tot
					end if
					
					arrayscal = maxloc(sum_prob_tot(1,:))
					indx_best = arrayscal(1)
					best_as = population(indx_best,:)
					my_results_population = population
					my_results_best = best_as
					my_results_overallsum = sum_prob_tot(1,:)
					opt_id = opt_idx_a(:,:,indx_best)
					
					do j = 1,expected_clusters
						optpar(j,:,:) = k_new(opt_id(j,:),:)
					end do
					
					my_results_parameter = optpar
					my_results_runtime_eachloop = enditeration-startiteration
					my_results_param_idx = opt_id
					
!				else
					
!					my_results_population = 0
!					opt_id = opt_idx_a(:,:,1)
					
!					do j = 1,expected_clusters
!						optpar(j,:,:) = k_new(opt_id(j,:),:)
!					end do
					
!					my_results_parameter = optpar
!					my_results_param_idx = opt_id
							
!				end if
				!my_results_population = population
				!print*, my_results_population
				!print*, "END Fortran"
				
				!print*, my_distance(1,:)
				!print*, clusterprobabilities(1,:)						
				 
				 !print*, my_results_runtime_eachloop
				!open(unit = 50, file = "/home/R/Fortran files/Greedy/my_results_population.txt")
				!write(50,*) my_results_population
				
				!open(unit = 60, file = "/home/R/Fortran files/Greedy/my_results_runtime_eachloop.txt")
				!write(60,*) my_results_runtime_eachloop
				
!----------------------------------------------------------------------------------------------!
!----------------------------Second Greedy run (optimize = 0)----------------------------------!
!----------------------------------------------------------------------------------------------!
				
		else
				
		do jj = 1,loops
		call cpu_time(startiteration(jj))
		!jj = 1
			as(:) = as_all(jj,:)
			
			
			max_mRNA_counts = 200
			num_cells_mRNA(:,:) = 0
			my_break = 0
				do iterations = 1,max_iter
		!print*, "Welcome3 to Fortran!" 
					call mysort(as,nvars,sorted_as,sorted_idx,as_sorted)
				!print*, "Welcome4 to Fortran!" 	
					do i = 1,nvars
					data_sorted(i,:) = mRNA_all(sorted_idx(i),:)
					end do
					
					do i = 1,expected_clusters
					clusters(i) = i
					end do
					
					j = 0
					
					do i = 1,nvars
						if (sorted_as(i)<sorted_as(i+1)) then
						j = j+1
						last_idx(j) = i
						end if
						last_idx(expected_clusters) = nvars
					end do
					
					num_clusters = size(clusters)
					
					bounds(:,2) = last_idx
					bounds(1,1) = 1
					
					do i = 2,num_clusters
						bounds(i,1) = bounds(i-1,2) + 1
					end do
					
					X = transpose(log_p_mat_ma)
					
					opt_idx_clusters(:,:) = 0
					logP_each_gene_each_clust(:,:) = 0
					
					do clust = 1,num_clusters
				
					IF( ALLOCATED(cells_in_each_cluster) )  DEALLOCATE( cells_in_each_cluster ) 						
					allocate(cells_in_each_cluster((bounds(clust,2)-bounds(clust,1)+1),n_genes))
				
					cells_in_each_cluster = data_sorted(bounds(clust,1):bounds(clust,2),:)	
					
					if (bounds(clust,1) == bounds(clust,2)) then
						print*, "Error: single cell in one cluster..."
					else
					!print*, "Welcome5 to Fortran!" 
						call myhist(expected_clusters,n_genes,num_cells_mRNA,max_mRNA_counts,bounds,clust,cells_in_each_cluster)
					
					end if		
					
					!print*, "Y before"
					 Y = num_cells_mRNA
					 !print*, "Y after"
					 Y(:,:) = Y(:,:) / sum(Y(:,1))
					 
					 !call multmat(X,Y,Z, n_genes, max_mRNA_counts)
					 !print*, "Welcome5.1 to Fortran!" 
					 Z = matmul(X,Y)
					 !print*, "Welcome5.2 to Fortran!"
					 do i = 1,n_genes
						Z_temp = Z(:,i)
						max_L(i) = maxval(Z_temp)
						!arrayscal is defined to save the array output, we then convert rank 1 array to scalar...
						arrayscal = maxloc(Z_temp)
						idx_max_L(i) = arrayscal(1)
					 end do
					 !print*, "Welcome5.3 to Fortran!"
					 opt_idx_clusters(clust,:) = idx_max_L
					 logP_each_gene_each_clust(clust,:) = max_L
					 
					 end do
									 
					 if (my_break == 1) exit
					 
					 idx_max_cell_prob(:) = 0
					 cell_prob(:,:) = 0.0
					 my_distance = cell_prob	
					 
					 do i = 1,nvars
						do clust = 1,num_clusters
							do j = 1,n_genes
								opt_param_each_gene = log_p_mat_ma(:,opt_idx_clusters(clust,j))
								cell_prob(i,clust) = cell_prob(i,clust) + opt_param_each_gene(mRNA_all(i,j)+1)
							end do
						end do
						!print*, "Welcome5.4 to Fortran!"
						sum_prob = maxval(cell_prob(i,:))
						arrayscal = maxloc(cell_prob(i,:))
						idx_max_cell_prob(i) = arrayscal(1)
						sum_prob_tot(1,jj) = sum_prob_tot(1,jj) + sum_prob
						maxrel = 1
						
					 end do
					 
					 do i = 1,nvars
					 sum_prob_max(i) = maxval(cell_prob(i,:))
					 arrayscal = maxloc(cell_prob(i,:))
					 idx_max_cell_prob(i) = arrayscal(1)
					 end do
					 
					 sum_prob_tot(1,jj) = sum(sum_prob_max)
					 
					 maxrel = 1
					 
							 
					 do i = 1,nvars
						if (as(i)==idx_max_cell_prob(i)) then
						comp(i) = 1
						else 
						comp(i) = 0
						end if
					 end do
						
					rel = sum(comp(:))
					
					!print*, "Welcome6 to Fortran!" 
					
					if (rel/nvars >= maxrel .or. optimize == 0) then
					
						population(jj,:) = as
						my_break = 1
						!print*, as
						!print*, shape(as)
						!if (optimize == 0) then
								do i = 1,expected_clusters
							!i = 1
									q = 0
									do j = 1,nvars

										if (as(j) == i) then
											q = q + 1
										end if
									end do
									!print*, q
									
									IF( ALLOCATED( aaa ) )  DEALLOCATE( aaa ) 	
									allocate(aaa(q))
									IF( ALLOCATED( bbb ) )  DEALLOCATE( bbb ) 	
									allocate(bbb(q,expected_clusters))
									aaa(:) = 0
									bbb(:,:) = 0
									call myfind(i,q,as,nvars,aaa)
									do j = 1,q
										bbb(j,:) = cell_prob(aaa(j),i)
									end do
									my_distance(aaa,:) = abs(bbb-cell_prob(aaa,:))
									
								end do
								log_clusterprobabilities = (log10(2.0))*cell_prob
								clusterprobabilities = 10**log_clusterprobabilities
								call myshape2(clusterprobabilities,expected_clusters,nvars,ccc)	
								clusterprobabilities = clusterprobabilities/ccc
								
								my_results_clusterprobabilities(:,:,jj) = clusterprobabilities
								my_results_distance(:,:,jj) = my_distance
								my_results_cell_prob(:,:,jj) = cell_prob
								my_results_best = as
							
						 !end if
					
					end if
					
					as = idx_max_cell_prob
					
				end do	 
					
				opt_idx_a(:,:,jj) = opt_idx_clusters
				
				call cpu_time(enditeration(jj))
				
				
				!print*, cell_prob(aaa,:)
				
				!do i = 1,q
				!print*, bbb(i,:)
				!end do
				
				!print*, jj
								
				end do
				
				
				
!				if (optimize==1) then
				
!					if (maxval(population(1,:)) /= maxval(population(2,:))) then
!						my_results_sum_prob_tot = sum_prob_tot
!					end if
					
!					arrayscal = maxloc(sum_prob_tot(1,:))
!					indx_best = arrayscal(1)
!					best_as = population(indx_best,:)
!					my_results_population = population
!					my_results_best = best_as
!					my_results_overallsum = sum_prob_tot(1,:)
!					opt_id = opt_idx_a(:,:,indx_best)
					
!					do j = 1,expected_clusters
!						optpar(j,:,:) = k_new(opt_id(j,:),:)
!					end do
					
!					my_results_parameter = optpar
!					my_results_runtime_eachloop = enditeration-startiteration
!					my_results_param_idx = opt_id
					
!				else
					
					my_results_population = 0
					opt_id = opt_idx_a(:,:,1)
					
					do j = 1,expected_clusters
						optpar(j,:,:) = k_new(opt_id(j,:),:)
					end do
					
					my_results_parameter = optpar
					my_results_param_idx = opt_id
							
				!end if
				
	
				!print*, as
				
				
				end if
				 	
end subroutine greedy_cabsel_f




!----------------------------------------------------------------------------------------------!
!----------------------------Sorting and tracking----------------------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine mysort(as,nvars,sorted_as,sorted_idx,as_sorted)

	integer i,j,temp(2),as(nvars), as_sorted(nvars,2), idx(nvars), sorted_idx(nvars), sorted_as(nvars)
	
	do i = 1,nvars
		idx(i) = i
	end do
	
	as_sorted(:,1) = as
	as_sorted(:,2) = idx
	
	do i = 1,nvars
		do j = i+1,nvars
			if (as_sorted(i,1) > as_sorted(j,1)) then
				temp = as_sorted(i,:)
				as_sorted(i,:) = as_sorted(j,:)
				as_sorted(j,:) = temp
			end if
		end do
	end do
	
	sorted_as(:) = as_sorted(:,1)
	sorted_idx(:) = as_sorted(:,2)
	
end subroutine mysort


!----------------------------------------------------------------------------------------------!
!--------------------------- # of cells with mRNA counts---------------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine myhist(expected_clusters,n_genes,num_cells_mRNA,max_mRNA_counts,bounds,clust,cells_in_each_cluster)

integer q,i,j,k,n_genes, clust, max_mRNA_counts, expected_clusters
integer num_cells_mRNA(201,n_genes), bounds(expected_clusters,2)
integer cells_in_each_cluster((bounds(clust,2)-bounds(clust,1)+1),n_genes),w
!print*, "myhist in"
!print*, clust
	do k = 1,n_genes
!	print*, "in k"
		do j = 0,max_mRNA_counts
!		print*, "in k & j"
		q = 0
			do i = 1,(bounds(clust,2)-bounds(clust,1)+1)
			!print*, "in k, k & i"
				if (cells_in_each_cluster(i,k) == j) then
				q = q + 1
				end if
				num_cells_mRNA(j+1,k) = q
			end do
		end do
	end do
	
!	print*, "myhist out"
!	print*, clust

end subroutine myhist

!----------------------------------------------------------------------------------------------!
!---------------------------------X and Y multiplication---------------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine multmat(X,Y,Z, n_genes, max_mRNA_counts)

integer i, j, k, n_genes, max_mRNA_counts

real*8 X(945,201), Y(201,n_genes), Z(945,n_genes), q

Z(:,:) = 0

do i = 1,945
	
	do k = 1,n_genes
	q = 0
		do j = 1,201
			
		q = q + X(i,j) * Y(j,k)
		
		end do
	Z(i,k) = q	
	end do
	
end do

end subroutine multmat

!----------------------------------------------------------------------------------------------!
!----------------------------------------------find--------------------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine myfind(i,q,as,nvars,aaa)

integer i,j,k,nvars,as(nvars),q,aaa(q)

k = 0

do j = 1,nvars
	if (as(j) == i) then
		k = k + 1
		aaa(k) = j
	end if
end do

end subroutine myfind		


!----------------------------------------------------------------------------------------------!
!-------------------------------------------Repmat2--------------------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine myshape2(clusterprobabilities,expected_clusters,nvars,ccc)

integer j, nvars, expected_clusters

real*8 clusterprobabilities(nvars,expected_clusters), ccc(nvars,expected_clusters)

do j = 1,nvars
ccc(j,:) = (sum(clusterprobabilities(j,:)))
end do

end subroutine myshape2
