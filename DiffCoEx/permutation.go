package main

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"strings"
	"time"

	"gonum.org/v1/gonum/mat"
)

// Helper function to calculate Spearman correlation
func spearmanCorrelation(x, y []float64) float64 {
	if len(x) != len(y) {
		panic("Input slices must have equal length")
	}

	n := len(x)

	// Create rank pairs
	type pair struct {
		value float64
		rank  float64
	}

	// Calculate ranks for x
	xPairs := make([]pair, n)
	for i, v := range x {
		xPairs[i] = pair{v, float64(i)}
	}
	sort.Slice(xPairs, func(i, j int) bool {
		return xPairs[i].value < xPairs[j].value
	})
	for i := range xPairs {
		xPairs[i].rank = float64(i + 1)
	}

	// Calculate ranks for y
	yPairs := make([]pair, n)
	for i, v := range y {
		yPairs[i] = pair{v, float64(i)}
	}
	sort.Slice(yPairs, func(i, j int) bool {
		return yPairs[i].value < yPairs[j].value
	})
	for i := range yPairs {
		yPairs[i].rank = float64(i + 1)
	}

	// Restore original order and calculate correlation
	xRanks := make([]float64, n)
	yRanks := make([]float64, n)
	for i := range xPairs {
		idx := int(xPairs[i].rank) - 1
		xRanks[idx] = float64(i + 1)
	}
	for i := range yPairs {
		idx := int(yPairs[i].rank) - 1
		yRanks[idx] = float64(i + 1)
	}

	// Calculate mean of ranks
	var sumX, sumY float64
	for i := 0; i < n; i++ {
		sumX += xRanks[i]
		sumY += yRanks[i]
	}
	meanX := sumX / float64(n)
	meanY := sumY / float64(n)

	// Calculate correlation
	var numerator, denomX, denomY float64
	for i := 0; i < n; i++ {
		dx := xRanks[i] - meanX
		dy := yRanks[i] - meanY
		numerator += dx * dy
		denomX += dx * dx
		denomY += dy * dy
	}

	return numerator / math.Sqrt(denomX*denomY)
}

// Helper function to convert mat.Dense column to []float64
func getColumn(m *mat.Dense, col int) []float64 {
	rows, _ := m.Dims()
	column := make([]float64, rows)
	mat.Col(column, col, m)
	return column
}

// Helper function to get columns matching a specific color
func getColumnsForColor(data *mat.Dense, colorh1C1C2 map[string]string, color string) [][]float64 {
	rows, cols := data.Dims()
	var columns [][]float64

	for i := 0; i < cols; i++ {
		if val, exists := colorh1C1C2[fmt.Sprintf("gene%d", i+1)]; exists && val == color {
			column := make([]float64, rows)
			mat.Col(column, i, data)
			columns = append(columns, column)
		}
	}
	return columns
}

// dispersionModule2Module calculates the dispersion value between two modules
func dispersionModule2Module(c1, c2 string, datC1, datC2 *mat.Dense, colorh1C1C2 map[string]string) float64 {
	if c1 == c2 {
		columnsC1 := getColumnsForColor(datC1, colorh1C1C2, c1)
		columnsC2 := getColumnsForColor(datC2, colorh1C1C2, c1)

		n := len(columnsC1)
		if n == 0 {
			return 0.0
		}

		var sumDifCorSquared float64
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				corC1 := spearmanCorrelation(columnsC1[i], columnsC1[j])
				corC2 := spearmanCorrelation(columnsC2[i], columnsC2[j])
				difCor := corC1 - corC2
				sumDifCorSquared += difCor * difCor
			}
		}

		denominator := float64(n*n-n) / 2.0
		return math.Sqrt((1.0 / denominator) * (sumDifCorSquared / 2.0))

	} else {
		columnsC1_1 := getColumnsForColor(datC1, colorh1C1C2, c1)
		columnsC1_2 := getColumnsForColor(datC1, colorh1C1C2, c2)

		columnsC2_1 := getColumnsForColor(datC2, colorh1C1C2, c1)
		columnsC2_2 := getColumnsForColor(datC2, colorh1C1C2, c2)

		n1 := len(columnsC1_1)
		n2 := len(columnsC1_2)
		if n1 == 0 || n2 == 0 {
			return 0.0
		}

		var sumDifCorSquared float64
		for i := 0; i < n1; i++ {
			for j := 0; j < n2; j++ {
				corC1 := spearmanCorrelation(columnsC1_1[i], columnsC1_2[j])
				corC2 := spearmanCorrelation(columnsC2_1[i], columnsC2_2[j])
				difCor := corC1 - corC2
				sumDifCorSquared += difCor * difCor
			}
		}

		return math.Sqrt((1.0 / float64(n1*n2)) * sumDifCorSquared)
	}
}

// generatePermutations creates a set of permuted indices
func generatePermutations(datC1, datC2 *mat.Dense, numPermutations int) [][]int {
	rows1, _ := datC1.Dims()
	rows2, _ := datC2.Dims()
	totalSamples := rows1 + rows2
	sampleSize := rows1

	permutations := make([][]int, numPermutations)
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	for i := 0; i < numPermutations; i++ {
		indices := make([]int, totalSamples)
		for j := range indices {
			indices[j] = j
		}

		for j := totalSamples - 1; j > 0; j-- {
			k := r.Intn(j + 1)
			indices[j], indices[k] = indices[k], indices[j]
		}

		permutations[i] = make([]int, sampleSize)
		copy(permutations[i], indices[:sampleSize])
	}

	return permutations
}

// scaleData scales the input data to have mean 0 and variance 1
func scaleData(data *mat.Dense) *mat.Dense {
	rows, cols := data.Dims()
	if rows == 0 {
		return data
	}

	scaledData := mat.NewDense(rows, cols, nil)

	for j := 0; j < cols; j++ {
		column := getColumn(data, j)

		// Calculate mean
		var sum float64
		for _, val := range column {
			sum += val
		}
		mean := sum / float64(rows)

		// Calculate variance
		var sumSquaredDiff float64
		for _, val := range column {
			diff := val - mean
			sumSquaredDiff += diff * diff
		}
		stdDev := math.Sqrt(sumSquaredDiff / float64(rows-1))

		// Scale the column
		for i := 0; i < rows; i++ {
			if stdDev != 0 {
				scaledData.Set(i, j, (data.At(i, j)-mean)/stdDev)
			} else {
				scaledData.Set(i, j, data.At(i, j)-mean)
			}
		}
	}

	return scaledData
}

// combineAndScaleData combines and scales datC1 and datC2
func combineAndScaleData(datC1, datC2 *mat.Dense) *mat.Dense {
	rows1, cols1 := datC1.Dims()
	rows2, cols2 := datC2.Dims()
	if cols1 != cols2 {
		panic("Matrices must have the same number of columns")
	}

	totalRows := rows1 + rows2
	combinedData := mat.NewDense(totalRows, cols1, nil)

	// Copy datC1
	for i := 0; i < rows1; i++ {
		for j := 0; j < cols1; j++ {
			combinedData.Set(i, j, datC1.At(i, j))
		}
	}

	// Copy datC2
	for i := 0; i < rows2; i++ {
		for j := 0; j < cols1; j++ {
			combinedData.Set(i+rows1, j, datC2.At(i, j))
		}
	}

	return scaleData(combinedData)
}

// permutationProcedureModule2Module calculates dispersion values using permuted data
func permutationProcedureModule2Module(permutation []int, d *mat.Dense, c1, c2 string, colorh1C1C2 map[string]string) float64 {
	rows, cols := d.Dims()

	// Create d1 from permuted indices
	d1 := mat.NewDense(len(permutation), cols, nil)
	for i, idx := range permutation {
		for j := 0; j < cols; j++ {
			d1.Set(i, j, d.At(idx, j))
		}
	}

	// Create d2 from remaining indices
	usedIndices := make(map[int]bool)
	for _, idx := range permutation {
		usedIndices[idx] = true
	}

	d2 := mat.NewDense(rows-len(permutation), cols, nil)
	j := 0
	for i := 0; i < rows; i++ {
		if !usedIndices[i] {
			for k := 0; k < cols; k++ {
				d2.Set(j, k, d.At(i, k))
			}
			j++
		}
	}

	return dispersionModule2Module(c1, c2, d1, d2, colorh1C1C2)
}

// Function to read the file and return a map
func readGeneColorFile(filename string) (map[string]string, error) {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Initialize the map to store gene-color pairs
	geneColorMap := make(map[string]string)

	// Create a scanner to read through the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		// Split each line by space, assuming that the first column is the gene and the second is the color
		parts := strings.Fields(line)
		if len(parts) == 2 {
			gene := parts[0]
			color := parts[1]
			// Add to the map
			geneColorMap[gene] = color
		} else {
			// Handle unexpected line format
			fmt.Printf("Skipping invalid line: %s\n", line)
		}
	}

	// Check for errors during scanning
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return geneColorMap, nil
}

/*

package main

import (
	"fmt"
)

func main() {
  // Example data as mat.Dense matrices
  datC1 := mat.NewDense(2, 2, []float64{1, 2, 3, 4})
  datC2 := mat.NewDense(2, 2, []float64{5, 6, 7, 8})



  // Specify the file name (change this to the path of your file)
	colorhC1C2 := "geneModuleColors.txt"

	// Read the file and get the gene-color map
	geneColorMap, err := readGeneColorFile(filename)
	if err != nil {
		fmt.Println("Error reading file:", err)
		return
	}

	// Output the map for verification
	fmt.Println("Gene-Color Map:")
	for gene, color := range geneColorMap {
		fmt.Printf("%s: %s\n", gene, color)
	}

    // Generate permutations
    numPermutations := 1000
    permutations := generatePermutations(datC1, datC2, numPermutations)

    // Scale and combine the data
    d := combineAndScaleData(datC1, datC2)

    // Compute dispersion matrix and null distribution
    uniqueColors := []string{"red", "blue"}
    dispersionMatrix := make([][]float64, len(uniqueColors))
    nullDistrib := make(map[string]map[string][]float64)

    for i, c1 := range uniqueColors {
        dispersionMatrix[i] = make([]float64, len(uniqueColors))
        nullDistrib[c1] = make(map[string][]float64)
        for j, c2 := range uniqueColors {
            dispersionMatrix[i][j] = dispersionModule2Module(c1, c2, datC1, datC2, colorh1C1C2)

            nullDistrib[c1][c2] = make([]float64, numPermutations)
            for k, perm := range permutations {
                nullDistrib[c1][c2][k] = permutationProcedureModule2Module(perm, d, c1, c2, colorh1C1C2)
            }
        }
    }

    // Create permutation summary
    permutationSummary := make([][]int, len(uniqueColors))
    for i := range permutationSummary {
        permutationSummary[i] = make([]int, len(uniqueColors))
        for j := range permutationSummary[i] {
            count := 0
            for _, val := range nullDistrib[uniqueColors[i]][uniqueColors[j]] {
                if val >= dispersionMatrix[i][j] {
                    count++
                }
            }
            permutationSummary[i][j] = count
        }
    }

    fmt.Println("Dispersion Matrix:", dispersionMatrix)
    fmt.Println("Permutation Summary:", permutationSummary)
}

*/
