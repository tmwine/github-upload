function lh_mx_out = Iso_LH_Graph_Chk(k_sym,node_mx,start_node)

%{

	given a graph and a starting node, this creates a canonical form for the graph, based on a "left-hand rule" exploration of all the nodes

	inputs:
	k_sym--number of symbols
	node_mx--this is the graph, in cylinder matrix form; a temporal graph matrix modded by the min_orbit length; rows of type [time now, node # at this time step, edge #, exit node (original num), this node (original num), time of exit node, num of exit node at that time]
	start_node--the node to start the left-hand walk on, the node in [time,node # at that time] format

	output:
	lh_mx_out--a matrix coding the graph's canonical form; it has k_sym+1 rows and (total # of nodes) columns

%}


%***
% setup
%***

node_org_num = start_node;
edge_stack = [];
conn_mx = []; % this has to start empty; its number of columns equals the number of new nodes discovered so far

node_num = 0; % which node we're on
% ***edges have form [this node # (new numbering), symbol/edge number [1-k_sym], exit node # (original numbering)]***
new_edge = [node_num,0,node_org_num]; % initializes routine; the source node # is null / "0", and symbol edge number is null / "0"
last_node = 0;

all_done = false;

edges_traversed = 0;



%***
% loop to perform the left-hand walk
%***

	% as the "left hand rule" walk is followed, the connection matrix, conn_mx, will be formed; it will have k_sym+1 rows and (total # of nodes) columns; the first column corresponds to the first node discovered, and so on; each row of a given column contains, in symbol # order (symbols are usually numbered from 1, so symbol 1 corresponds to 1st row, symbol 2 to 2nd, etc.) the nodes (if any) an edge w/ that symbol # leads to ("null" can be maybe -1, or 0), and the last row holds the original symbol numbers from the input graph (ie from node_mx); eg on k_sym=4, the columnar vector
	%	[10,-1,3,2,13]'
	%	indicates the node whose number is equal to that column's number in the matrix has exit edges a,c,d, going to nodes 10,3,2, respectively, and the node # is 13 (in node_mx)
while !all_done

	% does new_edge lead to a new node, or a node that's already been encountered?
	exit_node_abc = false; % this will be true if the exit node this edge leads to has already been covered
	if !isempty(conn_mx)
		nd_to = find(conn_mx(k_sym+1,:)==new_edge(3)); % ie looks amongst the original node numbers in the (k_sym+1)th row of conn_mx for a match--the index (column #) of this match is the node number (in new node numbering) of the match
		if !isempty(nd_to)
			exit_node_abc = true;
		endif
	endif

	if new_edge(1)>0 % count edges covered; don't count the first "null" edge leading to the first node to check
		edges_traversed++;
	endif

	if !exit_node_abc % ie if the node hasn't been encountered <-> if it's a new node

		node_num++;

		% update the connection matrix
		conn_mx(1:k_sym+1,node_num) = [-1*ones(k_sym,1);new_edge(3)]; % initialize this (node_num)th column with "nulls," and the (k_sym+1)th row with the (new) node number, in original numbering
		if last_node != 0
			conn_mx(new_edge(2),last_node) = node_num; % update the connection matrix node that lead to this one
		endif

		% a new node has been discovered; store its topmost (lowest symbol value) exit edge, and store the other edges in the edge stack (if there are others)
		% edges are stored as row vectors: [this node # (new numbering), symbol/edge number [1-k_sym], exit node # (original numbering)]
		ix_exit_edges = find(node_mx(:,5)==new_edge(3)); % originating node number (under original numbering) is in column 5 of node_mx
		tmp_len = length(ix_exit_edges);
		if tmp_len>1 % the new node has more than one exit edge; the left hand walk only follows the "topmost" edge (lowest-numbered symbol); we'll push (store) the other edges onto the stack, to return to (pull from) when the route encounters a node that's already been encountered
			tmp_block_mx = [node_num*ones(tmp_len,1),node_mx(ix_exit_edges,3),node_mx(ix_exit_edges,4)]; % each row has this node's number (the one w/ exit edges being listed) (in the new numbering), the symbol (edge) number, and the original node number the edge leads to
			[dum,srt_ix] = sort(tmp_block_mx(:,2),"descend"); % find row indices of matrix sorted by symbol #s (descending order, starting from top)
			ee_block_mx = tmp_block_mx(srt_ix,:);
			new_edge = ee_block_mx(end,:);
			ee_block_mx(end,:) = []; % pop this edge from the temporary edges matrix (edges as rows)
			edge_stack = [edge_stack; ee_block_mx]; % push the other edges onto the stack
		else
			new_edge = [node_num,node_mx(ix_exit_edges,3),node_mx(ix_exit_edges,4)];
		endif

		last_node = node_num;

	else % ie if the node has already been encountered / covered

		% record this edge in the connection matrix
		conn_mx(new_edge(2),last_node) = nd_to;

		if !isempty(edge_stack)
			% "backtrack"--pop next most "lefthanded" edge/path from the stack (the bottom row of edge_stack matrix)
			new_edge = edge_stack(end,:);
			edge_stack(end,:) = []; % remove this edge from the stack
			last_node = new_edge(1);
		else % here if this edge leads to a node that's already been checked, and there are no more edges to check / backtrack to
			all_done = true;
		endif

	endif

endwhile


if edges_traversed!=size(node_mx,1)
	% this is likely an error in the code
	printf("problem in left hand rule graph rendering\n");
endif

lh_mx_out = conn_mx;


endfunction
