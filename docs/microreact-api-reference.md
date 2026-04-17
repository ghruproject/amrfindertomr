# Microreact .microreact JSON Format Reference

Minimum working format for creating Microreact projects via API or file upload. Based on reverse-engineering working projects — the official docs don't cover this fully.

## API Endpoint

```
POST https://microreact.org/api/projects/create
Content-Type: application/json; charset=utf-8
Access-Token: <your-token>
```

Get your token from https://microreact.org/my-account/settings

## Minimum Working JSON

Tested and confirmed working as of April 2026. Every field below is required unless marked optional.

```json
{
  "schema": "https://microreact.org/schema/v1.json",
  "charts": {},
  "datasets": {
    "dataset-1": {
      "id": "dataset-1",
      "file": "data-file-1",
      "idFieldName": "id"
    }
  },
  "files": {
    "data-file-1": {
      "id": "data-file-1",
      "format": "text/csv",
      "name": "metadata.csv",
      "type": "data",
      "size": 95,
      "blob": "data:application/octet-stream;base64,<BASE64_ENCODED_CSV>"
    },
    "tree-file-1": {
      "id": "tree-file-1",
      "format": "text/x-nh",
      "name": "tree.nwk",
      "type": "tree",
      "size": 49,
      "blob": "data:application/octet-stream;base64,<BASE64_ENCODED_NEWICK>"
    }
  },
  "filters": {
    "paneFilters": [],
    "dataFilters": [],
    "chartFilters": [],
    "searchOperator": "includes",
    "searchValue": "",
    "selection": [],
    "selectionBreakdownField": null
  },
  "maps": {},
  "matrices": {},
  "meta": {
    "name": "My Project Name"
  },
  "networks": {},
  "notes": {},
  "panes": { "...see Panes section below..." },
  "slicers": {},
  "styles": {
    "coloursField": null,
    "colourPalettes": [],
    "defaultColour": "transparent",
    "defaultShape": "circle",
    "colourSettings": {},
    "labelsField": null,
    "legendDirection": "row",
    "shapesField": null,
    "shapePalettes": []
  },
  "tables": {
    "table-1": {
      "displayMode": "cosy",
      "hideUnselected": false,
      "dataset": "dataset-1",
      "file": "data-file-1",
      "title": "Metadata",
      "paneId": "table-1",
      "columns": [
        { "field": "id", "fixed": false },
        { "field": "column2", "fixed": false }
      ]
    }
  },
  "timelines": {},
  "trees": {
    "tree-1": { "...see Trees section below..." }
  },
  "views": []
}
```

## Key Gotchas

1. **All top-level keys are required** — even empty ones (`charts: {}`, `maps: {}`, etc.). Missing keys cause "editundefined" errors.

2. **`idFieldName` must match a CSV column** — this is the column that links metadata rows to tree tip labels.

3. **Files use `blob`, not `url`** — embed data as `data:application/octet-stream;base64,<data>`. Using `url` doesn't work for API uploads.

4. **Files need `type` and `size`** — `type` is `"data"` for CSV or `"tree"` for Newick. `size` is byte length.

5. **Table needs `columns`, `paneId`, `file`, and `dataset`** — without `columns`, the table panel renders empty. Without `paneId`, the panel shows "editundefined".

6. **Pane tab IDs must be unique** — duplicate IDs (e.g. two tabs with `id: "tree-1"`) crash the client with "each node must have a unique id".

7. **Pane tab IDs must match defined panels** — a tab with `id: "tree-1"` requires a matching `trees.tree-1` definition. A tab with `id: "table-1"` requires `tables.table-1`.

8. **`styles` and `filters` need full structure** — empty `{}` for these causes rendering issues. Include all subfields even if null/empty.

## Files

```json
{
  "id": "data-file-1",
  "format": "text/csv",
  "name": "metadata.csv",
  "type": "data",
  "size": 1234,
  "blob": "data:application/octet-stream;base64,aWQsZ2VuZTEKc2FtcGxlMSx5ZXMK"
}
```

| Field | Values | Notes |
|-------|--------|-------|
| `format` | `text/csv`, `text/x-nh` | CSV for data, x-nh for Newick |
| `type` | `"data"`, `"tree"` | Must match the file content |
| `size` | integer | Byte length of the original (unencoded) content |
| `blob` | `data:application/octet-stream;base64,...` | Base64-encoded file content |

## Datasets

Links a data file to the project and specifies the ID column.

```json
{
  "id": "dataset-1",
  "file": "data-file-1",
  "idFieldName": "id"
}
```

- `file` must reference a file ID from the `files` section
- `idFieldName` must exactly match a column header in the CSV

## Tables

```json
{
  "displayMode": "cosy",
  "hideUnselected": false,
  "dataset": "dataset-1",
  "file": "data-file-1",
  "title": "Metadata",
  "paneId": "table-1",
  "columns": [
    { "field": "id", "fixed": false },
    { "field": "gene1", "fixed": false }
  ]
}
```

- `columns` must list every CSV column you want displayed
- `paneId` must match the tab `id` in the panes layout
- Both `dataset` and `file` references are needed

## Trees

Full tree configuration (all fields needed for proper rendering):

```json
{
  "alignLabels": true,
  "blockHeaderFontSize": 13,
  "blockPadding": 0,
  "blocks": [],
  "blockSize": 14,
  "branchLengthsDigits": 4,
  "controls": false,
  "fontSize": 16,
  "hideOrphanDataRows": false,
  "ids": null,
  "internalLabelsFilterRange": [0, 100],
  "internalLabelsFontSize": 13,
  "lasso": false,
  "nodeSize": 14,
  "path": null,
  "roundBranchLengths": true,
  "scaleLineAlpha": true,
  "showBlockHeaders": true,
  "showBlockLabels": false,
  "showBranchLengths": false,
  "showEdges": true,
  "showInternalLabels": false,
  "showLabels": true,
  "showLeafLabels": false,
  "showPiecharts": true,
  "showShapeBorders": true,
  "showShapes": true,
  "styleLeafLabels": false,
  "styleNodeEdges": false,
  "subtreeIds": null,
  "type": "rc",
  "title": "Tree",
  "labelField": "id",
  "file": "tree-file-1"
}
```

- `type`: `"rc"` = rectangular, `"rd"` = radial, `"dg"` = diagonal
- `labelField`: which CSV column to use for tip labels
- `file`: must reference a tree file ID

## Panes Layout

Uses [FlexLayout](https://github.com/nickhudkins/FlexLayout) model. The layout defines how panels are arranged.

### Side-by-side (tree left, table right)

```json
{
  "model": {
    "global": {
      "splitterSize": 2,
      "tabEnableClose": false,
      "tabSetHeaderHeight": 1,
      "tabSetTabStripHeight": 1,
      "tabSetMinWidth": 160,
      "tabSetMinHeight": 160,
      "borderMinSize": 160,
      "borderBarSize": 20,
      "borderEnableDrop": false
    },
    "borders": [
      {
        "type": "border",
        "size": 240,
        "location": "right",
        "children": [
          { "type": "tab", "id": "--mr-legend-pane", "name": "Legend", "component": "Legend", "enableClose": false, "enableDrag": false },
          { "type": "tab", "id": "--mr-selection-pane", "name": "Selection", "component": "Selection", "enableClose": false, "enableDrag": false },
          { "type": "tab", "id": "--mr-history-pane", "name": "History", "component": "History", "enableClose": false, "enableDrag": false },
          { "type": "tab", "id": "--mr-views-pane", "name": "Views", "component": "Views", "enableClose": false, "enableDrag": false }
        ]
      }
    ],
    "layout": {
      "type": "row",
      "children": [
        {
          "type": "tabset",
          "weight": 60,
          "children": [
            { "type": "tab", "id": "tree-1", "name": "Tree", "component": "Tree" }
          ]
        },
        {
          "type": "tabset",
          "weight": 40,
          "children": [
            { "type": "tab", "id": "table-1", "name": "Metadata", "component": "Table" }
          ]
        }
      ]
    }
  }
}
```

### Stacked (tree on top, table below)

Wrap each tabset in a `"type": "row"` node — the outer layout type stays `"row"` but the nested rows stack vertically:

```json
"layout": {
  "type": "row",
  "children": [
    {
      "type": "row",
      "children": [
        { "type": "tabset", "weight": 60, "children": [
          { "type": "tab", "id": "tree-1", "name": "Tree", "component": "Tree" }
        ]}
      ]
    },
    {
      "type": "row",
      "children": [
        { "type": "tabset", "weight": 40, "children": [
          { "type": "tab", "id": "table-1", "name": "Metadata", "component": "Table" }
        ]}
      ]
    }
  ]
}
```

### Available components

| Component | Panel type |
|-----------|-----------|
| `Tree` | Phylogenetic tree |
| `Table` | Metadata table |
| `Map` | Geographic map |
| `Timeline` | Timeline |
| `Matrix` | Distance matrix |
| `Chart` | Chart |
| `Legend` | Legend (border only) |
| `Selection` | Selection (border only) |
| `History` | History (border only) |
| `Views` | Views (border only) |

## API Response

```json
{
  "id": "unique-project-id",
  "url": "https://microreact.org/project/unique-project-id"
}
```

## Python Example

```python
import base64, json, requests

csv_data = open("metadata.csv").read()
tree_data = open("tree.nwk").read()

csv_b64 = base64.b64encode(csv_data.encode()).decode()
tree_b64 = base64.b64encode(tree_data.encode()).decode()

project = {
    "schema": "https://microreact.org/schema/v1.json",
    # ... full JSON as above, with blob fields set to:
    # f"data:application/octet-stream;base64,{csv_b64}"
    # f"data:application/octet-stream;base64,{tree_b64}"
}

resp = requests.post(
    "https://microreact.org/api/projects/create",
    json=project,
    headers={
        "Content-Type": "application/json; charset=utf-8",
        "Access-Token": "your-token-here",
    },
)
print(resp.json()["url"])
```
