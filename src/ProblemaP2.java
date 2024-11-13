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
			llenarRed(n, d, celulas, peptidos, red);
			int [] nodoMetodo2=edmonsKarpNodoMaximoFlujo(red,celulas);
			int celulaBloqueada=nodoMetodo2[0];
			int flujoMax=nodoMetodo2[1];
			int flujoMin=nodoMetodo2[2];
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
			for (DirectedWeightedEdge edge: residual.getAdjacents(u)) {
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
			for (DirectedWeightedEdge edge: residual.getAdjacents(u)) {
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
			for(DirectedWeightedEdge edge : red.getAdjacents(nodo)) {
				int dest = edge.getDest();
				residual.addEdge(nodo,dest,edge.getWeight());
				residual.addEdge(dest, nodo, 0);
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
				for (DirectedWeightedEdge edge: residual.getAdjacents(u)) {
					if(edge.getDest() == v) {
						capacity = edge.getWeight();
						flow = edge.getFlow();					}
				}
				pathFlow = Math.min(pathFlow, capacity-flow);
				v=u;
			}
			
			v = sumidero;
			while(v!=source) {
				int u = parent[v];
				
				for (DirectedWeightedEdge edge: residual.getAdjacents(u)) {
					if(edge.getDest()==v) {
						edge.setFlow(edge.getFlow()+pathFlow);
					}
				}
				
				for (DirectedWeightedEdge edge: residual.getAdjacents(v)) {
					if(edge.getDest()==u) {
						edge.setFlow(edge.getFlow()-pathFlow);
					}
				}
				
				v=u;
				
			}
			
			flujo += pathFlow;
			
			HashMap<Integer,ArrayList<DirectedWeightedEdge>> crossing = new HashMap<Integer,ArrayList<DirectedWeightedEdge>>();
			HashMap<Integer, Integer> crossingFlows = new HashMap<Integer, Integer>();
        	int anterior=-909;

			boolean visitados[] = bfsResidual(residual, source, sumidero, parent);
			for (int o: residual.getSet()) {
				for (DirectedWeightedEdge p: residual.getAdjacents(o)) {
					if (visitados[o] && !visitados[p.getDest()]) {
						crossing.computeIfAbsent(o, k -> new ArrayList<>()).add(p);
				        int newFlow = crossingFlows.getOrDefault(0, 0) + p.getFlow();
				        crossingFlows.put(o, newFlow);
				        if (crossingFlows.getOrDefault(o,Integer.MIN_VALUE)>=maxCrossFlow) {
				        	if (celulas.containsKey(p.getDest())){
					        	maxCrossFlow = crossingFlows.getOrDefault(o,Integer.MIN_VALUE);
					        	int tipo = celulas.get(p.getDest()).getTipo();
					        	anterior=parent[p.getDest()];
					        	maxNode = p.getDest();
					        	if (tipo!=2) maxNode=anterior;
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
				
				//double manhattanDist = Math.sqrt((celulaOrigenX-celulaDestinoX)^2 + (celulaOrigenY-celulaDestinoY)^2);
				int manhattanDist = Math.abs(celulaOrigenX-celulaDestinoX)+Math.abs(celulaOrigenY-celulaDestinoY);
				if (manhattanDist<=d && manhattanDist!=0) {
				
					if ((celulaOrigenTipo==1 && celulaDestinoTipo==2) || (celulaOrigenTipo==2 && celulaDestinoTipo==3)){
					   
					    int mensajes = calcularMensajes(celulas.get(i), celulas.get(j), peptidos);
						red.addEdge(celulaOrigenId, celulaDestinoId, mensajes);
					}
					
					if (celulaOrigenTipo==2 && celulaDestinoTipo==2){
						//verifico que no exista
						boolean exists = false;
					    for (DirectedWeightedEdge edge : red.getAdjacents(celulaDestinoId)) {
					        if (edge.getDest() == celulaOrigenId) {
					            exists = true;
					            break;
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
		private Map<Integer, List<DirectedWeightedEdge>> adjList;
		public DirectedWeightedGraph() {
			adjList = new HashMap<>();
		}
		
		public void addEdge(int source, int dest, int weight) {
		    adjList.putIfAbsent(source, new ArrayList<>());
		    adjList.get(source).add(new DirectedWeightedEdge(dest, weight));
	    }
		
		
		public List<DirectedWeightedEdge> getAdjacents(int source){
			return adjList.getOrDefault(source, new ArrayList<>());
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
	


