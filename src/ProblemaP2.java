import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.ArrayList;
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
			//System.out.println("voy a procesar n "+n+" d: "+d);

			//celula[][] malla = new celula[n][n];
			int peptideCount = 0;

			for (int c=0;c<n;c++) {
				
				linea = lector.readLine();
				
				String[] campos = linea.split("\\s+");
				
				int id = Integer.parseInt(campos[0]);
				int x = Integer.parseInt(campos[1]);
				//System.out.println("estoy leyendo id "+id+" x: "+x);
				int y = Integer.parseInt(campos[2]);
				int tipo = Integer.parseInt(campos[3]);
				HashSet<Integer> peptidosCelula = new HashSet<>();
				
				celula celula = instance.new celula(id,x,y,tipo,peptidosCelula);
				celulas.put(id, celula);
				//malla[celula.getX()-1][celula.getY()-1]=celula;
				
				for (int e=4;e<campos.length;e++) {
					if(!peptidos.containsKey(campos[e])) {
						peptidos.put(campos[e], peptideCount);
						peptidosCelula.add(peptidos.get(campos[e]));
						peptideCount++;
						//System.out.println("celula "+celula.getId()+" peptido: "+campos[e]+" codigo: "+peptidos.get(campos[e]));
					}	
					else {
						peptidosCelula.add(peptidos.get(campos[e]));
						//System.out.println("celula "+celula.getId()+" peptido: "+campos[e]+" codigo: "+peptidos.get(campos[e]));
					}
				 }
			}
			llenarRed(n, d, celulas, peptidos, red);
			int flujoMax = edmonsKarp(red);
			int flujoMin = Integer.MAX_VALUE;
			int celulaBloqueada = -1;
			//System.out.println(celulas.keySet());
			for (int key:celulas.keySet()) {
				if (celulas.get(key).getTipo()!=2) continue;
				int celulaId = celulas.get(key).getId();
				int flujoBloqueo = edmonsKarpBloqueo(red,celulaId);
				//System.out.println("bloquee: "+celulaId+" flujo: "+flujoBloqueo);
				if (flujoBloqueo<=flujoMin) {
					flujoMin = flujoBloqueo;
					celulaBloqueada = celulaId;
				}
			}
			System.out.println(Integer.toString(celulaBloqueada)+" "+Integer.toString(flujoMax)+" "+Integer.toString(flujoMin));
		}
	}
	
	private static boolean bfs(DirectedWeightedGraph residual, int source, int sumidero, int[] parent) {
		//System.out.println("mi tamaño del residual es: "+sizeresidual);
		parent[source]=-1;
		boolean[] visitados = new boolean[sumidero+1];
		visitados[source]=true;
		//HashMap<Integer, Boolean> visitados = new HashMap<Integer, Boolean>();
		Queue<Integer> cola = new LinkedList<>();
		cola.add(source);
		
		//for (int i: residual.getSet()) {
		//    visitados[i]=false;
		//}
		
		while(!cola.isEmpty()) {
			int u = cola.poll();
			for (DirectedWeightedEdge edge: residual.getAdjacents(u)) {
				int v = edge.getDest();
				//System.out.println("bfs ya lo visite: "+visitados.getOrDefault(v,false));
				//System.out.println("bfs mi peso: "+edge.getWeight());
				//System.out.println("bfs mi flow: "+edge.getFlow());

				if(!visitados[v] && (edge.getWeight() - edge.getFlow() > 0)) {
					parent[v]=u;
					//parent.put(v, u);
					//visitados.put(v, true);
					visitados[v]=true;
					cola.add(v);
					//System.out.println("bfs v: "+v+" bfs sumidero: "+sumidero);
					if (v == sumidero) {
						//System.out.println("bfs encontro true");
						return true;
					}
				}
			}
		}
		return false;
	}
	
	private static DirectedWeightedGraph crearResidualBloqueo(DirectedWeightedGraph red, int bloqueado) {
		DirectedWeightedGraph residual = new DirectedWeightedGraph();
		for (int source: red.getSet()) {
			if (source==bloqueado) continue;
			for(DirectedWeightedEdge edge : red.getAdjacents(source)) {
				int dest = edge.getDest();
				if(dest==bloqueado) continue;
				residual.addEdge(source,dest,edge.getWeight());
				residual.addEdge(dest, source, 0);
			}
		}
		return residual;
	}
	
	private static int edmonsKarpBloqueo(DirectedWeightedGraph red, int bloqueado) {
		int source = 0;
		int sumidero = red.getSize();
		//System.out.println("bloqueado: "+bloqueado);
		int flujo = 0;
		DirectedWeightedGraph residual = crearResidualBloqueo(red, bloqueado);
		int[] parent = new int[red.getSize()+1];
		//HashMap<Integer,Integer>parent = new HashMap<Integer,Integer>();
		
		while(bfs(residual,source,sumidero,parent)) {
			int pathFlow = Integer.MAX_VALUE;
			int v = sumidero;
			
			while(v != source) {
				int u = parent[v];
				int capacity = Integer.MAX_VALUE;
				int flow = 0;
				for (DirectedWeightedEdge edge: residual.getAdjacents(u)) {
					if(edge.getDest()==v) {
						capacity = edge.getWeight();
						flow = edge.getFlow();
					}
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
		}	
 		return flujo;
	}
	
	
	private static DirectedWeightedGraph crearResidual(DirectedWeightedGraph red) {
		DirectedWeightedGraph residual = new DirectedWeightedGraph();
		for (int nodo: red.getSet()) {
			for(DirectedWeightedEdge edge : red.getAdjacents(nodo)) {
				int dest = edge.getDest();
				//System.out.println("destino: "+ dest+ " peso: "+edge.getWeight());
				residual.addEdge(nodo,dest,edge.getWeight());
				residual.addEdge(dest, nodo, 0);
			}
		}
		return residual;
	}
	
	private static int edmonsKarp(DirectedWeightedGraph red) {
		int source = 0;
		int sumidero = red.getSize();
		//System.out.println("sumidero sin bloqueo: "+sumidero);
		int flujo = 0;
		DirectedWeightedGraph residual = crearResidual(red);
		//System.out.println("size de grafo residual en iteración :"+residual.getSize()+" size grafo: "+red.getSize());
		int[] parent = new int[red.getSize()+1];
		//HashMap<Integer,Integer> parent = new HashMap<Integer,Integer>();

		while(bfs(residual, source, sumidero, parent)) {
			//System.out.println("entre a while ek");
			int pathFlow = Integer.MAX_VALUE;
			int v = sumidero;
			//System.out.println("v en edmons karp interación: "+Integer.toString(v));
			while(v != source) {
				//System.out.println("v en iteracion ek: "+ v);

				int u = parent[v];
				int capacity = 0;
				int flow = 0;
				for (DirectedWeightedEdge edge: residual.getAdjacents(u)) {
					//System.out.println("edge en iteracion ek: "+ edge.getDest());
					if(edge.getDest() == v) {
						capacity = edge.getWeight();
						flow = edge.getFlow();
						//System.out.println("flow en iteración ek: "+flow);
					}
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
		}	
 		return flujo;
	}
	
	
	private static void llenarRed(int n, int d, HashMap<Integer, celula> celulas, HashMap<String, Integer> peptidos, DirectedWeightedGraph red) {
		int adicion = 1;
		//System.out.println(celulas.keySet());
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
				
				int manhattanDist = Math.abs(celulaOrigenX-celulaDestinoX)+Math.abs(celulaOrigenY-celulaDestinoY);
				if (manhattanDist<=d && manhattanDist!=0) {
				
					if ((celulaOrigenTipo==1 && celulaDestinoTipo==2) || (celulaOrigenTipo==2 && celulaDestinoTipo==3)){
					   
					    int mensajes = calcularMensajes(celulas.get(i), celulas.get(j), peptidos);
						red.addEdge(celulaOrigenId, celulaDestinoId, mensajes);
					    //System.out.println("creo conexion: "+celulaOrigenId+" "+celulaDestinoId+ " peso: "+mensajes);
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
							
							//añado eje de ida
							red.addEdge(celulaOrigenId, celulaDestinoId, mensajes);
						    //System.out.println("creo conexion: "+celulaOrigenId+" "+celulaDestinoId+ " peso: "+mensajes);
			            	//añado eje de vuelta con nodo auxiliar
							//String auxiliary = Integer.toString(celulaOrigenId) + "0" + Integer.toString(celulaDestinoId);
			            	int auxiliaryId = n+adicion;
			            	red.addEdge(celulaDestinoId,auxiliaryId,mensajes);
			            	red.addEdge(auxiliaryId, celulaOrigenId, mensajes);
						    //System.out.println("creo conexion: "+celulaDestinoId+" "+celulaOrigenId+ " peso: "+mensajes);

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
			    //System.out.println("creo conexion: "+0+" "+celulaId+ " peso: "+Integer.MAX_VALUE);
				//System.out.println("adyacentes a las primeras: "+red.getAdjacents(-1).size());
				//System.out.println("tamaño red final: "+Integer.toString(red.getSize()));
			}
		}
		
		int t=n+adicion;
		
		for (int k:celulas.keySet()) {
			int celulaTipo = celulas.get(k).getTipo();
			int celulaId = celulas.get(k).getId();
			if (celulaTipo==3) {
				red.addEdge(celulaId,t,Integer.MAX_VALUE);
			    //System.out.println("creo conexion: "+celulaId+" "+t+ " peso: "+Integer.MAX_VALUE);

			}
		}
	}
				
	private static int calcularMensajes(celula origen, celula destino, HashMap<String, Integer> peptidos) {
		int peso = 0;
		for (Integer peptidoOrigen : origen.getPeptidos()) {
			if (destino.getPeptidos().contains(peptidoOrigen)) {
				//System.out.println("coinciden peptidos: "+peptidoOrigen+ "mensajes: "+peso);
				peso++;
			}
		}
		//System.out.println("mensajes: "+peso);
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
		    //depende de que me responda Tomás sobre la entrada

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
	


