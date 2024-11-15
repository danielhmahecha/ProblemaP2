import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.lang.Math;

public class ProblemaP2 {
	
	public static void main(String[] args) throws Exception{
		ProblemaP2 instance = new ProblemaP2();
		
		BufferedReader lector = new BufferedReader(new InputStreamReader(System.in));
		int casos = Integer.parseInt(lector.readLine());
		
		for (int p=0;p<casos;p++) {
			String linea = lector.readLine();

			HashMap<String, Integer> peptidos = new HashMap<String,Integer>();			
			HashMap<Integer, celula> celulas = new HashMap<Integer, celula>();
			DirectedWeightedGraph red = new DirectedWeightedGraph();
			
			String[] items = linea.split("\\s+");
			
			int n = Integer.parseInt(items[0]);
			int d = Integer.parseInt(items[1]);

			int peptideCount = 0;

			for (int c=0;c<n;c++) {
				
				linea = lector.readLine();
				
				String[] campos = linea.split("\\s+");
				
				int id = Integer.parseInt(campos[0]);
				int x = Integer.parseInt(campos[1]);
				int y = Integer.parseInt(campos[2]);
				int tipo = Integer.parseInt(campos[3]);
				HashSet<Integer> peptidosCelula = new HashSet<>();
				
				celula celula = instance.new celula(id,x,y,tipo,peptidosCelula);
				celulas.put(id, celula);
				
				for (int e=4;e<campos.length;e++) {
					if(!peptidos.containsKey(campos[e])) {
						peptidos.put(campos[e], peptideCount);
						peptidosCelula.add(peptidos.get(campos[e]));
						peptideCount++;
					}	
					else {
						peptidosCelula.add(peptidos.get(campos[e]));
					}
				 }
			}
			if (d/n > 0.5) d= (int) Math.round(n*0.7);
			llenarRed(n, d, celulas, peptidos, red);
			int [] nodoMetodo2=edmonsKarpNodoMaximoFlujo(red,celulas);
			int celulaBloqueada=nodoMetodo2[0];
			int flujoMax=nodoMetodo2[1];
			int flujoMin=nodoMetodo2[2];
			//int flujoPush=getPushRelabel(red);
			//int flujoPush=0;
			System.out.println(Integer.toString(celulaBloqueada)+" "+Integer.toString(flujoMax)+" "+Integer.toString(flujoMin));
		}
	}
	
	private static boolean bfs(DirectedWeightedGraph residual, int source, int sumidero, int[] parent) {
		parent[source]=-1;
		boolean[] visitados = new boolean[sumidero+1];
		visitados[source]=true;
		Queue<Integer> cola = new LinkedList<>();
		cola.add(source);

		while(!cola.isEmpty()) {
			int u = cola.poll();
			for (DirectedWeightedEdge edge: residual.getAdjacents(u).values()) {
				int v = edge.getDest();
	
				if(!visitados[v] && (edge.getWeight() - edge.getFlow() > 0)) {
					parent[v]=u;
					visitados[v]=true;
					cola.add(v);
					if (v == sumidero) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
	private static boolean[] bfsResidual(DirectedWeightedGraph residual, int source, int sumidero, int[] parent) {
		parent[source]=-1;
		boolean[] visitados = new boolean[sumidero+1];
		Arrays.fill(visitados, false);
		visitados[source]=true;
		Queue<Integer> cola = new LinkedList<>();
		cola.add(source);
	
		while(!cola.isEmpty()) {
			int u = cola.poll();
			for (DirectedWeightedEdge edge: residual.getAdjacents(u).values()) {
				int v = edge.getDest();

				if(!visitados[v] && (edge.getWeight() - edge.getFlow() > 0)) {
					parent[v]=u;

					visitados[v]=true;
					cola.add(v);
					if (v == sumidero) {
						return visitados;
					}
				}
			}
		}
		return visitados;
	}
		
	private static DirectedWeightedGraph crearResidual(DirectedWeightedGraph red) {
		DirectedWeightedGraph residual = new DirectedWeightedGraph();
		for (int nodo: red.getSet()) {
			for(DirectedWeightedEdge edge : red.getAdjacents(nodo).values()) {
				int dest = edge.getDest();
				residual.addEdge(nodo,dest,edge.getWeight());
				residual.addEdge(dest, nodo, 0);
				//System.out.println("añadí "+nodo+" "+dest+" "+edge.getWeight());
				//System.out.println("añadí "+dest+" "+nodo+" "+0);

			}
		}
		return residual;
	}
	
	private static DirectedWeightedGraph crearResidualPush(DirectedWeightedGraph red) {
		DirectedWeightedGraph residual = new DirectedWeightedGraph();
		for (int nodo: red.getSet()) {
			for(DirectedWeightedEdge edge : red.getAdjacents(nodo).values()) {
				int dest = edge.getDest();
				residual.addEdge(nodo,dest,edge.getWeight());
				//System.out.println("añadí "+nodo+" "+dest+" "+edge.getWeight());
				//System.out.println("añadí "+dest+" "+nodo+" "+0);

			}
		}
		return residual;
	}
	private static int[] edmonsKarpNodoMaximoFlujo(DirectedWeightedGraph red, HashMap<Integer, celula> celulas ) {
		int[] answer = new int[3];
		int maxNode = 0;
		int maxCrossFlow = Integer.MIN_VALUE;
		int source = 0;
		int sumidero = red.getSize();
		int flujo = 0;
		DirectedWeightedGraph residual = crearResidual(red);
		int[] parent = new int[red.getSize()+1];

		while(bfs(residual, source, sumidero, parent)) {
			int pathFlow = Integer.MAX_VALUE;
			int v = sumidero;
			while(v != source) {
				
				int u = parent[v];
				int capacity = 0;
				int flow = 0;
				DirectedWeightedEdge edge = residual.getEdge(u, v);
				capacity = edge.getWeight();
				flow = edge.getFlow();					
				pathFlow = Math.min(pathFlow, capacity-flow);
				v=u;
			}
			
			v = sumidero;
			while(v!=source) {
				int u = parent[v];
				DirectedWeightedEdge edge = residual.getEdge(u,v);
				edge.setFlow(edge.getFlow()+pathFlow);
				
				DirectedWeightedEdge vedge = residual.getEdge(v, u);
				vedge.setFlow(vedge.getFlow()-pathFlow);
				v=u;
			}
			
			flujo += pathFlow;
			
			HashMap<Integer,ArrayList<DirectedWeightedEdge>> crossing = new HashMap<Integer,ArrayList<DirectedWeightedEdge>>();
			HashMap<Integer, Integer> crossingFlows = new HashMap<Integer, Integer>();
        	int anterior=-909;

			boolean visitados[] = bfsResidual(residual, source, sumidero, parent);
			for (int o: residual.getSet()) {
				for (int p: residual.getSet()) {
					if (residual.getEdge(o,p)!=null && (visitados[o] && !visitados[p])) {
						crossing.computeIfAbsent(o, k -> new ArrayList<>()).add(residual.getEdge(o, p));
				        int newFlow = crossingFlows.getOrDefault(0, 0) + residual.getEdge(o, p).getFlow();
				        crossingFlows.put(o, newFlow);
				        if (crossingFlows.getOrDefault(o,Integer.MIN_VALUE)>=maxCrossFlow) {
				        	if (celulas.containsKey(residual.getEdge(o, p).getDest())){
					        	maxCrossFlow = crossingFlows.getOrDefault(o,Integer.MIN_VALUE);
					        	int tipo = celulas.get(residual.getEdge(o, p).getDest()).getTipo();
					        	if (tipo!=2) maxNode=parent[residual.getEdge(o, p).getDest()];
					        	else maxNode = residual.getEdge(o, p).getDest();

				        	}
				        }
					}
				}
			}
		}	
 		answer[1]=flujo;
		answer[2]=flujo-maxCrossFlow;
 		answer[0]=maxNode;
		return answer;
	}
	
	private static void llenarRed(int n, int d, HashMap<Integer, celula> celulas, HashMap<String, Integer> peptidos, DirectedWeightedGraph red) {
		int adicion = 1;
		for (int i:celulas.keySet()) {
			for (int j:celulas.keySet()) {
				
				int celulaOrigenX = celulas.get(i).getX();
				int celulaOrigenY = celulas.get(i).getY();
				
				int celulaDestinoX = celulas.get(j).getX();
				int celulaDestinoY = celulas.get(j).getY();
				
				int celulaOrigenTipo = celulas.get(i).getTipo();
				int celulaDestinoTipo = celulas.get(j).getTipo();
				
				int celulaOrigenId = celulas.get(i).getId();
				int celulaDestinoId = celulas.get(j).getId();
				
				double euclidean = Math.sqrt((celulaOrigenX-celulaDestinoX)*(celulaOrigenX-celulaDestinoX) + (celulaOrigenY-celulaDestinoY)*(celulaOrigenY-celulaDestinoY));
				//System.out.println("euclidean dist "+manhattanDist);
				//int manhattanDist = Math.abs(celulaOrigenX-celulaDestinoX)+Math.abs(celulaOrigenY-celulaDestinoY);
				//System.out.println("euclidean dist "+manhattanDist+" manhattan "+manhattan);

				if (euclidean<=d && euclidean!=0) {
				
					if ((celulaOrigenTipo==1 && celulaDestinoTipo==2) || (celulaOrigenTipo==2 && celulaDestinoTipo==3)){
					   
					    int mensajes = calcularMensajes(celulas.get(i), celulas.get(j), peptidos);
						red.addEdge(celulaOrigenId, celulaDestinoId, mensajes);
					}
					
					if (celulaOrigenTipo==2 && celulaDestinoTipo==2){
						//verifico que no exista
						boolean exists = false;
						if (red.containsKey(celulaDestinoId)) {
						    for (DirectedWeightedEdge edge : red.getAdjacents(celulaDestinoId).values()) {
						        if (edge.getDest() == celulaOrigenId) {
						            exists = true;
						            break;
						        }
						    }
						}
					    if (!exists) {
							int mensajes = calcularMensajes(celulas.get(i), celulas.get(j), peptidos);
							
							red.addEdge(celulaOrigenId, celulaDestinoId, mensajes);
			            	int auxiliaryId = n+adicion;
			            	red.addEdge(celulaDestinoId,auxiliaryId,mensajes);
			            	red.addEdge(auxiliaryId, celulaOrigenId, mensajes);

			            	adicion++;
					    }
					}	
				}
			}
		}
		
		for (int k:celulas.keySet()) {
			int celulaTipo = celulas.get(k).getTipo();
			int celulaId = celulas.get(k).getId();
			//Creo super source
			if (celulaTipo==1) {
				red.addEdge(0, celulaId, Integer.MAX_VALUE);
			}
		}
		
		int t=n+adicion;
		
		for (int k:celulas.keySet()) {
			int celulaTipo = celulas.get(k).getTipo();
			int celulaId = celulas.get(k).getId();
			if (celulaTipo==3) {
				red.addEdge(celulaId,t,Integer.MAX_VALUE);

			}
		}
	}
				
	private static int calcularMensajes(celula origen, celula destino, HashMap<String, Integer> peptidos) {
		int peso = 0;
		for (Integer peptidoOrigen : origen.getPeptidos()) {
			if (destino.getPeptidos().contains(peptidoOrigen)) {
				peso++;
			}
		}
		return peso;
	}
	
	/*
	private static boolean push (int node,DirectedWeightedGraph residual,int[] altura, int[] exceso, Queue<Integer> activeNodes) {
		System.out.println("entré a push");
		for (DirectedWeightedEdge e: residual.getAdjacents(node).values()) {
			System.out.println("Altura node "+altura[node]+" altura dest "+altura[e.getDest()]);
			if((altura [node] > altura[e.getDest()]) && (e.getFlow()!=e.getWeight())){
				int flow = Math.min(e.getWeight()-e.getFlow(), exceso[node]);
				exceso[node] -= flow;
				exceso[e.getDest()] += flow;
				e.setFlow(e.getFlow()+flow);
				System.out.println("nuevo exceso: "+e.getDest()+" "+exceso[e.getDest()]);
				boolean a = true;
				
				if (residual.containsKey(e.getDest())) {
					if (residual.getAdjacents(e.getDest()).containsKey(node)) {			
						residual.getEdge(e.getDest(), node).setFlow(-flow);
						a=false;
					} else {
						residual.addEdge(e.getDest(), node, 0);
						residual.getEdge(e.getDest(), node).setFlow(-flow);
						a=false;
					}
				}
				
				if (a) {
					residual.addEdge(e.getDest(), node, 0);
					residual.getEdge(e.getDest(), node).setFlow(-flow);
				}
				if (exceso[e.getDest()] > 0 && !activeNodes.contains(e.getDest())) activeNodes.add(e.getDest());

				return true;
				}
			}
		
		return false;
	}
	
	private static  void relabel (int node,DirectedWeightedGraph residual,int[] altura, int[] exceso) {
		int minAltura = Integer.MAX_VALUE;
		if(residual.containsKey(node)){
		for (DirectedWeightedEdge e: residual.getAdjacents(node).values()) {
			if((e.getFlow()!=e.getWeight()) && (altura[e.getDest()]<minAltura)) {
				minAltura=altura[e.getDest()];
				altura[node]=minAltura+1;
				System.out.println(" nueva altura "+altura[node]);
				if(exceso[e.getDest()] > 0 && e.getDest() != 0 && e.getDest() != residual.getSize()) {
					activeNodes.add(e.getDest());
					
				}
			}
		}
		
		
	}
	
	private static int getPushRelabel(DirectedWeightedGraph red) {
		//creo nuevo residual
		DirectedWeightedGraph residual = crearResidualPush(red);
		
		//inicializo altura y exceso
		int[] altura = new int[residual.getSize()+1];
		Arrays.fill(altura, 0);
		int[] exceso = new int[residual.getSize()+1];
		Arrays.fill(exceso,0);
		
		
		//lleno activeNodes
		Queue<Integer> activeNodes= new LinkedList<Integer>();
        
		//preflow
		altura[0] = red.getSize()+1;
		for (DirectedWeightedEdge e: residual.getAdjacents(0).values()) {
			e.setFlow(e.getWeight());
			System.out.println("mi flow es: "+e.getFlow());
			exceso[e.getDest()] += e.getFlow();
			residual.addEdge(e.getDest(), 0, 0);
			if (exceso[e.getDest()] > 0) activeNodes.add(e.getDest());
			
			}
        //ciclo

		while (!activeNodes.isEmpty()) {
			int activeNode = activeNodes.poll();
			System.out.println("entré a while nodo: "+activeNode);
			if(!push(activeNode, residual, altura, exceso, activeNodes)) {
				relabel(activeNode,residual,altura,exceso);
			}
		}
		
		int maxFlow = exceso[residual.getSize()];
		return maxFlow;
	}*/
	
	
	class celula{
		private int id;
		private int x;
        private int y;
        private int tipo;
        private HashSet<Integer> peptidos;

        public celula(int id, int x, int y, int tipo, HashSet<Integer> peptidos) {
            this.id = id;
        	this.x = x;
            this.y = y;
            this.tipo = tipo;
            this.peptidos = peptidos;
        }
        
        public int getId() {
        	return id;
        }
        
        public int getX() {
            return x;
        }

        public int getY() {
            return y;
        }

        public int getTipo() {
            return tipo;
        }

        public HashSet<Integer> getPeptidos() {
            return peptidos;
        }
	}
	
	public static class DirectedWeightedEdge {
		int dest;
		int weight;
		int flow;
		
		public DirectedWeightedEdge(int dest, int weight) {
			super();
			this.dest = dest;
			this.weight = weight;
			this.flow = 0;
		}
		public int getDest() {
			return dest;
		}
		
		public int getWeight() {
			return weight;
		}
		
		public int getFlow() {
			return flow;
		}
		public void setFlow(int flow){
			this.flow=flow;
		}
	}
	
	public static class DirectedWeightedGraph {
		private Map<Integer, Map<Integer, DirectedWeightedEdge>> adjList;
		public DirectedWeightedGraph() {
			adjList = new HashMap<>();
		}
		
		public void addEdge(int source, int dest, int weight) {
		    adjList.putIfAbsent(source, new HashMap<Integer, DirectedWeightedEdge>());
		    adjList.get(source).put(dest, new DirectedWeightedEdge(dest, weight));
	    }
		
		public DirectedWeightedEdge getEdge (int source, int dest) {
			return adjList.get(source).get(dest);
		}
		
		public Map<Integer, DirectedWeightedEdge> getAdjacents(int source){
			return adjList.get(source);
		}
		
		public boolean containsKey (int key) {
			return adjList.containsKey(key);
		}
		
		public int getSize() {
			int size = adjList.size();
			return size;
		}
		
		public Set<Integer> getSet() {
			Set<Integer> keySet = adjList.keySet();
			return keySet;
		}
	}
	}
	


