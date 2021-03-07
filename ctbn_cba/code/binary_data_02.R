n=2
name="binary_data_02"
vars3_data = list(vars=data.frame("Name"=c("X","Y","Z"),"Value"=rep(n,3),stringsAsFactors = FALSE),
		  n_iter=10,
		  edge_prob=0.2,
		  time_end=100,
		  nsample=300,
		  name=name)

vars4_data = list(vars=data.frame("Name"=c("X","Y","Z","Q"),"Value"=rep(n,4),stringsAsFactors = FALSE),
		  n_iter=10,
		  edge_prob=0.2,
		  time_end=100,
		  nsample=300,
		  name=name)

vars5_data = list(vars=data.frame("Name"=c("X","Y","Z","Q","V"),"Value"=rep(n,5),stringsAsFactors = FALSE),
		  n_iter=10,
		  edge_prob=0.2,
		  time_end=100,
		  nsample=300,
		  name=name)

vars6_data = list(vars=data.frame("Name"=c("X","Y","Z","Q","V","A"),"Value"=rep(n,6),stringsAsFactors = FALSE),
		  n_iter=10,
		  edge_prob=0.2,
		  time_end=100,
		  nsample=300,
		  name=name)

vars10_data = list(vars=data.frame("Name"=c("X","Y","Z","Q","V","A","B","C","D","E"),"Value"=rep(n,10),stringsAsFactors = FALSE),
		   n_iter=3,
		   edge_prob=0.2,
		   time_end=100,
		   nsample=300,
		   name=name)

vars15_data = list(vars=data.frame("Name"=c("X","Y","Z","Q","V","A","B","C","D","E","F","G","H","I","J"),"Value"=rep(n,15),stringsAsFactors = FALSE),
		   n_iter=3,
		   edge_prob=0.2,
		   time_end=100,
		   nsample=300,
		   name=name)
vars20_data = list(vars=data.frame("Name"=c("X","Y","Z","Q","V","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O"),"Value"=rep(n,20),stringsAsFactors = FALSE),
 		   n_iter=3,
		   edge_prob=0.2,
		   time_end=100,
		   nsample=300,
		   name=name)

iter_list = list(vars3_data,
		 vars4_data,
		 vars5_data,
		 vars6_data,
		 vars10_data,
		 vars15_data,
		 vars20_data
		 )

