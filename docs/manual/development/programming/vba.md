# Vba

**vba** or visual-basic for applications is a scripting/programming language and development environment embedded within the Office suite software from microsoft. The syntax is equivalent to visual-basic, but the libraries and functions give direct access to the objects of the Office environment.

## Examples

### Sum an array

*Pt()* is a VBA array

```
Dim Pt() as Double
ReDim Pt(1 to 10)
Cells(2, 2) = WorksheetFunction.Sum(Pt)
```

### Write an array in the worksheet

The quick way to write a VBA array in the excel worksheet is by using the *Range* method:

```
Dim Pt() as Double
ReDim Pt(1 to 10)
Range(Cells(1, 1), Cells(1, 10)) = Pt
```

## Sample code

### Delete text

Subroutine that deletes text in excel.

```
'Rotina que apaga todos os valores na tabela de Excel
Public Sub ApagaTexto()

    Range(Cells(1, 5), Range("E1:E1").End(xlToRight).End(xlDown)).ClearContents

End Sub
```

### Delete charts

Subroutine that deletes all charts present in excel sheet

```
'Rotina que apaga os gráficos existentes
Public Sub apagaGrafico()

    'Iterador de tipo "gráfico"
    Dim Chart As ChartObject

    'Para cada gráfico dentro da colecção de gráficos da ActiveSheet ...
    For Each Chart In ActiveSheet.ChartObjects

        '... apaga o gráfico.
        Chart.Delete

    Next

End Sub
```

### Create chart

Simple subroutine that creates a chart

```
'Rotina que faz gráfico
Sub fazGrafico(nRows As Long, nCols As Long)

    ActiveSheet.Shapes.AddChart.Select
    ActiveChart.SetSourceData Source:=Range(Cells(2, 6), Cells(nRows + 1, nCols + 5)), _
        PlotBy:=xlRows
    ActiveChart.ChartType = xlLine
    ActiveChart.SeriesCollection(1).XValues = Range(Cells(1, 6), Cells(1, nCols + 5))

End Sub
```

Complex subroutine that creates a nice chart

```
'Rotina que cria um grafico na mesma folha que a chama
Sub makegraphic(lineSpan, colSpan As Integer, _
                strtLine, strtCol As Integer, _
                yScale As Single)

    Dim graphic As Chart                'Cria uma variavel de grafico
    Dim i As Single

    'Cria um gráfico dentro da active sheet
    Set graphic = ActiveSheet.ChartObjects.Add _
        (Left:=60, Width:=500, Top:=100, Height:=300).Chart

    With graphic

        .ChartType = xlLine
        .SetSourceData _
            Source:=Range(Cells(strtLine, strtCol), Cells(strtLine - 1 + lineSpan, strtCol - 1 + colSpan)), _
            PlotBy:=xlRows

        'define as abcissas
        .SeriesCollection(1).XValues = Range(Cells(strtLine - 1, strtCol), Cells(strtLine - 1, strtCol - 1 + colSpan))

        'Põe os eixos
        .HasTitle = True
        .ChartTitle.Characters.Text = "Crescimento de flores numa estufa"
        .Legend.Clear

        With .Axes(xlCategory, xlPrimary)
            .HasTitle = True
            .AxisTitle.Characters.Text = "Domínio (m)"
        End With

        With .Axes(xlValue)
            .HasTitle = True
            .AxisTitle.Characters.Text = "Tamanho (cm)"
            .TickLabels.NumberFormat = "0.0"
            .MaximumScale = yScale
            .MinimumScale = 0#
        End With

    End With

End Sub
```

### Thomas algorithm

The Thomas algorithm efficiently solves linear systems of equations who possess a tridiagonal matricial representation.

```
'----------------------------------------------------------------------------------------------------
'Este algoritmo resolve a equação:
'concent (i+1)*e_coef(i) + concent(i) * p_coef(i) + concent(i-1)*w_coef(i) = Ti(i)
'Referência (wikipedia): http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
'----------------------------------------------------------------------------------------------------
Sub Thomas(ByRef e_coef() As Double, _
            ByRef p_coef() As Double, _
            ByRef w_coef() As Double, _
            ByRef TI() As Double)

Dim PF() As Double
Dim QF() As Double

Dim U As Long
Dim L As Long

Dim j As Long

U = UBound(e_coef)
L = LBound(e_coef)

ReDim PF(L To U) As Double
ReDim QF(L To U) As Double

 PF(L) = -e_coef(L) / p_coef(L)
 QF(L) = TI(L) / p_coef(L)

  For j = L + 1 To U
   PF(j) = -e_coef(j) / (p_coef(j) + w_coef(j) * PF(j - 1))
   QF(j) = (TI(j) - w_coef(j) * QF(j - 1)) / (p_coef(j) + w_coef(j) * PF(j - 1))
  Next

 TI(U) = QF(U)

 For j = U - 1 To L Step -1
  TI(j) = PF(j) * TI(j + 1) + QF(j)
 Next

End Sub
```

## External links

- [Class of environmental modeling given at IST, 2008-2009](https://fenix.ist.utl.pt/homepage/ist146730/modelacao-ambiental-2008-2009)
