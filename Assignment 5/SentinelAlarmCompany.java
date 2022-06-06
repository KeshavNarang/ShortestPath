import java.io.*;
import java.util.*;
public class SentinelAlarmCompany
{
	public static void main (String [] args) throws IOException
	{
		Scanner inputReader = new Scanner (new FileReader (new File ("input.txt")));
		
		int numberVertices = inputReader.nextInt();
		Graph allNodes = new Graph (numberVertices);
		Graph onlySeven = new Graph (numberVertices);
		Graph onlyFourteen = new Graph (numberVertices);
		Graph onlyTwentyFive = new Graph (numberVertices);
		Graph onlyTwentyEight = new Graph (numberVertices);
		Graph onlyForty = new Graph (numberVertices);
		
		while (inputReader.hasNextLine())
		{
			int start = inputReader.nextInt();
			int end = inputReader.nextInt();
			int weight = inputReader.nextInt();
			
			allNodes.addEdge(new Edge(start, end, weight));
			onlySeven.addEdge(new Edge(start, end, weight));
			onlyFourteen.addEdge(new Edge(start, end, weight));
			onlyTwentyFive.addEdge(new Edge(start, end, weight));
			onlyTwentyEight.addEdge(new Edge(start, end, weight));
			onlyForty.addEdge(new Edge(start, end, weight));
		}
		
		// Part 1 - All the stations (7, 14, 25, 28, and 40) are connected to an imaginary node (50) with zero weight

		System.out.println("When all five stations are connected: ");
		Dijkstra shortestPath1 = new Dijkstra(allNodes, 50);
		
		// Part 2 - Keep only one station, and run Dijkstra's. The other four (as well as the imaginary node) are removed.
		
		System.out.println("When seven is the only station: ");
		onlySeven.removeNode(14);
		onlySeven.removeNode(25);
		onlySeven.removeNode(28);
		onlySeven.removeNode(40);
		onlySeven.removeNode(50);
		Dijkstra shortestPath2 = new Dijkstra(onlySeven, 7);
		
		System.out.println("When fourteen is the only station: ");
		onlyFourteen.removeNode(7);
		onlyFourteen.removeNode(25);
		onlyFourteen.removeNode(28);
		onlyFourteen.removeNode(40);
		onlyFourteen.removeNode(50);
		Dijkstra shortestPath3 = new Dijkstra(onlyFourteen, 14);
		
		System.out.println("When twenty-five is the only station: ");
		onlyTwentyFive.removeNode(7);
		onlyTwentyFive.removeNode(14);
		onlyTwentyFive.removeNode(28);
		onlyTwentyFive.removeNode(40);
		onlyTwentyFive.removeNode(50);
		Dijkstra shortestPath4 = new Dijkstra(onlyTwentyFive, 25);
		
		System.out.println("When twenty-eight is the only station: ");
		onlyTwentyEight.removeNode(7);
		onlyTwentyEight.removeNode(14);
		onlyTwentyEight.removeNode(25);
		onlyTwentyEight.removeNode(40);
		onlyTwentyEight.removeNode(50);
		Dijkstra shortestPath5 = new Dijkstra(onlyTwentyEight, 28);
		
		System.out.println("When forty is the only station: ");
		onlyForty.removeNode(7);
		onlyForty.removeNode(15);
		onlyForty.removeNode(25);
		onlyForty.removeNode(28);
		onlyForty.removeNode(50);
		Dijkstra shortestPath6 = new Dijkstra(onlyForty, 40);
		
		System.out.println("So if only one station can be kept, it should be 7.");
	}
}

class Dijkstra
{
	private Edge [] edgeTo;
	private double [] distTo;
	private IndexMinPQ <Double> pq;
	public Dijkstra (Graph G, int s)
	{
		edgeTo = new Edge[G.numberVertices];
		distTo = new double[G.numberVertices];
		pq = new IndexMinPQ<Double>(G.numberVertices);
		
		for (int v = 0; v < G.numberVertices; v++)
		{
			distTo[v] = Double.POSITIVE_INFINITY;
		}
		
		distTo[s] = 0.0;
		pq.insert(s, 0.0);
		
		while (!pq.isEmpty())
		{
			relax(G, pq.delMin());
		}
		
		double maxDistance = 0;
		for (double distance : distTo)
		{
			if (distance == Double.POSITIVE_INFINITY)
			{
				continue;
			}
			
			maxDistance = Math.max(maxDistance, distance);
		}			
		System.out.println("The max response time to any customer is " + maxDistance + ". \n");
	}
	
	private void relax (Graph G, int v)
	{
		for (Edge e : G.edges.get(v))
		{
			int w = e.other(v);
			if (distTo[w] > distTo[v] + e.weight)
			{
				distTo[w] = distTo[v] + e.weight;
				edgeTo[w] = e;
				
				if (pq.contains(w))
				{
					pq.change(w, distTo[w]);
				}

				else 
				{
					pq.insert(w, distTo[w]);
				}
			}
		}
	}
	
	public double distTo(int v)
	{ 
		return distTo[v];
	}
	
	public boolean hasPathTo(int v)
	{ 
		return distTo[v] < Double.POSITIVE_INFINITY; 
	}
	
	public Iterable <Edge> pathTo(int v)
	{
		if (!hasPathTo(v))
		{			
			return null;
		}
		Stack <Edge> path = new Stack <Edge>();
		for (Edge e = edgeTo[v]; e != null; e = edgeTo[e.other(v)])
		{
			path.push(e);
		}
		return path;
	}
}

class Graph
{
	public int numberVertices;
	public int numberEdges;

	public ArrayList <ArrayList<Edge>> edges;

	public Graph (int numberVertices)
	{
		this.numberVertices = numberVertices;
		this.numberEdges = 0;

		edges = new ArrayList <ArrayList<Edge>> ();
		for (int index = 0; index < numberVertices; index++)
		{
			edges.add(new ArrayList <Edge> ());
		}
	}

	public void addEdge (Edge edge)
	{
		int one = edge.one;
		int two = edge.two;
		edges.get(one).add(edge);
		edges.get(two).add(edge);
		numberEdges++;
	}
	
	public void removeNode (int node) // If you want to remove seven from the graph
	{
		for (Edge connected : edges.get(node))
		{
			int other = connected.other(node); // Go to every node seven is directly linked to
			
			Edge find = null;
			
			for (Edge next : edges.get(other))
			{
				if (next.other(other) == node)
				{
					find = next;
					break;
				}
			}
			
			edges.get(other).remove(find);
		}
		edges.get(node).clear(); // Clear the adjacency list of seven
	}
}

class Edge implements Comparable <Edge>
{
	public int one;
	public int two;
	public int weight;

	public Edge (int one, int two, int weight)
	{
		this.one = one;
		this.two = two;
		this.weight = weight;
	}

	public int other (int vertex)
	{
		if (one == vertex)
		{
			return two;
		}
		else
		{
			return one;
		}
	}

	public int compareTo (Edge other)
	{
		if (this.weight < other.weight)
		{
			return -1;
		}
		else if (this.weight > other.weight)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	public String toString ()
	{
		return one + " - " + two;
	}
}


class IndexMinPQ <Key extends Comparable<Key>> implements Iterable <Integer> 
{
    private int maxN;        
    private int n;        
    private int[] pq;     
    private int[] qp;      
    private Key[] keys;     

    public IndexMinPQ (int maxN) 
	{
        if (maxN < 0)
		{
			throw new IllegalArgumentException();
		}
		
        this.maxN = maxN;
        n = 0;
        keys = (Key[]) new Comparable[maxN + 1]; 
        pq = new int[maxN + 1];
        qp = new int[maxN + 1];   
		
        for (int i = 0; i <= maxN; i++)
		{
            qp[i] = -1;
		}
    }

    public boolean isEmpty() 
	{
        return n == 0;
    }

    public boolean contains (int i) 
	{
        validateIndex(i);
        return qp[i] != -1;
    }

    public int size() 
	{
        return n;
    }

    public void insert (int i, Key key) 
	{
        validateIndex(i);
        if (contains(i)) 
		{
			throw new IllegalArgumentException("index is already in the priority queue");
		}
		
        n++;
        qp[i] = n;
        pq[n] = i;
        keys[i] = key;
        swim(n);
    }

    public int minIndex() 
	{
        if (n == 0) 
		{
			throw new NoSuchElementException("Priority queue underflow");
		}
		
        return pq[1];
    }

    public Key minKey() 
	{
        if (n == 0)
		{
			throw new NoSuchElementException("Priority queue underflow");
		}
        return keys[pq[1]];
    }

    public int delMin() 
	{
        if (n == 0) 
		{
			throw new NoSuchElementException("Priority queue underflow");
		}
		
        int min = pq[1];
        exch(1, n--);
        sink(1);
        assert min == pq[n+1];
        qp[min] = -1;        
        keys[min] = null;    
        pq[n+1] = -1;        
        return min;
    }

    public Key keyOf (int i) 
	{
        validateIndex(i);
        if (!contains(i)) 
		{
			throw new NoSuchElementException("index is not in the priority queue");
		}
        else return keys[i];
    }

    public void changeKey (int i, Key key) 
	{
        validateIndex(i);
        if (!contains(i)) 
		{
			throw new NoSuchElementException("index is not in the priority queue");
		}
        keys[i] = key;
        swim(qp[i]);
        sink(qp[i]);
    }


    @Deprecated
    public void change(int i, Key key) 
	{
        changeKey(i, key);
    }

    public void decreaseKey (int i, Key key) 
	{
        validateIndex(i);
        if (!contains(i)) 
		{
			throw new NoSuchElementException("index is not in the priority queue");
		}
        if (keys[i].compareTo(key) == 0)
		{
            throw new IllegalArgumentException("Calling decreaseKey() with a key equal to the key in the priority queue");
		}
        if (keys[i].compareTo(key) < 0)
		{
            throw new IllegalArgumentException("Calling decreaseKey() with a key strictly greater than the key in the priority queue");
		}
        keys[i] = key;
        swim(qp[i]);
    }

    public void increaseKey (int i, Key key) 
	{
        validateIndex(i);
        if (!contains(i)) 
		{
			throw new NoSuchElementException("index is not in the priority queue");
		}
        if (keys[i].compareTo(key) == 0)
		{
            throw new IllegalArgumentException("Calling increaseKey() with a key equal to the key in the priority queue");
		}
        if (keys[i].compareTo(key) > 0)
		{
            throw new IllegalArgumentException("Calling increaseKey() with a key strictly less than the key in the priority queue");
		}
        keys[i] = key;
        sink(qp[i]);
    }

    public void delete (int i) 
	{
        validateIndex(i);
        if (!contains(i)) 
		{
			throw new NoSuchElementException("index is not in the priority queue");
		}
		
        int index = qp[i];
        exch(index, n--);
        swim(index);
        sink(index);
        keys[i] = null;
        qp[i] = -1;
    }

    private void validateIndex (int i) 
	{
        if (i < 0)
		{
			throw new IllegalArgumentException("index is negative: " + i);
		}
		
        if (i >= maxN) 
		{	
			throw new IllegalArgumentException("index >= capacity: " + i);
		}
    }

    private boolean greater (int i, int j) 
	{
        return keys[pq[i]].compareTo(keys[pq[j]]) > 0;
    }

    private void exch (int i, int j) 
	{
        int swap = pq[i];
        pq[i] = pq[j];
        pq[j] = swap;
        qp[pq[i]] = i;
        qp[pq[j]] = j;
    }

    private void swim (int k) 
	{
        while ((k > 1) && greater(k/2, k)) 
		{
            exch(k, k/2);
            k = k/2;
        }
    }

    private void sink (int k) 
	{
        while (2*k <= n) 
		{
            int j = 2*k;
			
            if (j < n && greater(j, j+1)) 
			{
				j++;
			}
			
            if (!greater(k, j))
			{
				break;
			}
			
            exch(k, j);
            k = j;
        }
    }

    public Iterator<Integer> iterator() 
	{ 
		return new HeapIterator(); 
	}

    private class HeapIterator implements Iterator<Integer> 
	{
        private IndexMinPQ<Key> copy;
		
        public HeapIterator() 
		{
            copy = new IndexMinPQ <Key> (pq.length - 1);
            for (int i = 1; i <= n; i++)
			{
                copy.insert(pq[i], keys[pq[i]]);
			}
        }

        public boolean hasNext()  
		{ 
			return !copy.isEmpty();                    
		}
		
        public void remove()  
		{
			throw new UnsupportedOperationException(); 
		}

        public Integer next() 
		{
            if (!hasNext())
			{
				throw new NoSuchElementException();
			}
            return copy.delMin();
        }
    }
}
