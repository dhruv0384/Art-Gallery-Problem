# Computational Geometry Project

## Art Gallery Problem & Mobile Edge Guard Extension

This project implements the **Art Gallery Problem** using computational geometry techniques, and extends it to the **Mobile Guard Problem**, where guards can patrol polygon edges or diagonals instead of being static at vertices. The system performs polygon triangulation, 3-coloring, vertex guard placement, and an experimental mobile edge guard algorithm.  

Screenshots of example runs are available in the `sample_results/` folder.

---

## 🚀 How to Run

1. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the main script**
    ```bash
    python main.py
    ```
3. A GUI will prompt you for the number of polygon vertices (`n`).  
The program will then:

    - Generate a simple polygon  
    - Partition into monotone polygons  
    - Triangulate  
    - Apply 3-coloring  
    - Show minimum vertex guards  
    - Show mobile edge guards  

---

## ✨ Features

- Polygon generation & DCEL representation  
- Monotone partitioning and triangulation  
- 3-coloring of triangulated dual graph (for vertex guard solution)  
- Vertex guard placement (O’Rourke’s theorem, ⌊n/3⌋ bound)  
- Mobile edge guard extension (experimental algorithm using diagonals and edge coverage)  
- Step-by-step visualization in a single GUI frame sequence  

---

## 📘 References

- O’Rourke, J. (1987). *Art Gallery Theorems and Algorithms*.  
- Recent algorithmic results on mobile edge guards (see `report.pdf` for details).  
