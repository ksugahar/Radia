# Onshape assembly export without a personal API key

Use an existing fixed-version Assembly tab URL with the path
`/documents/<document>/v/<version>/e/<element>`. Preserve the selected
configuration when recording the export. Agents may view documents and export
to local storage. They must ask the user to create a missing fixed version.
Version, workspace, branch, tab and copy creation, edits, and exports stored
in the document are document writes.

## Export and inspect

For a mechanism, start with Simscape Multibody on MATLAB R2026a. The exporter
uses the user's Onshape login and OAuth authorization.

```matlab
xmlFile = smexportonshape(url, 'FolderPath', outputDir);
```

Use an existing writable output directory. Use forward slashes in Windows
paths passed through shell commands to avoid escaping errors. Allow around
40 minutes for a large export. The output is a Multibody XML file and STEP
files for the part geometries. Retain the XML with its STEP files.

Use `build123d_probe` with `query="step_readback"` for the declared length unit,
volume in mm³ and centroid in mm. A part STEP gives a centroid in its part
frame. Compose the XML instance transforms along the full hierarchy to
express that centroid in the top assembly frame. The XML files inspected here
declare translation in metres. Convert it to millimetres before adding it to
the rotated centroid. The comparison used row-major 3 × 3 rotations. The
simple reference assembly's nonidentity rotation is symmetric, so that case
cannot distinguish row-major and column-major parsing.

For the geometry-only route, download a single assembly STEP from the same
fixed version and configuration in the Onshape UI. Preserve the assembly
axes and check the expected instance count. The UI export excludes hidden
instances unless **Include hidden instances** is enabled. Record that setting
and compare it with the visibility of each occurrence in the fixed version.
Keep the version URL, selected
configuration, export settings and file SHA-256 with the result.

If the workflow needs structured part or configuration metadata from the
REST API, use the Onshape adapter as a second choice. Users supply
`ONSHAPE_ACCESS_KEY` and `ONSHAPE_SECRET_KEY` in the server environment.

## Observations from 2026-11-23 and 2026-11-24

| Case | Observed result |
| --- | --- |
| Base-and-arm fixed version | Two unique part STEP files and three placed instances. Part-frame STEP centroids matched the XML within 4 × 10⁻⁷ mm. XML mass divided by STEP volume implied 2700 kg/m³ for both parts. The UI assembly STEP independently reproduced the three placed volumes and centroids to numerical precision. |
| Two-part nested assembly | The XML retained the parent and child structure. The fixed version and Main placed a part differently even though its volume was the same. Pin the version when comparing placements. |
| Configurable Assembly V28 | The default configuration exported 27 instances, an XML file and part STEP files. Screw STEP files with heights of 14, 10 and 7 mm shared a Part UID in the XML. All 12 screw instances referenced that UID. The initial UI STEP matched 11 screws. One screw appeared as hidden in the instance list. With **Include hidden instances** enabled, the bulk STEP matched all 12 screw placements one-to-one by volume, overall length and XML-transformed centroid. The maximum centroid distance was 2.56 × 10⁻⁷ mm. Only the download option changed. No visibility or suppression edits were needed. |

A STEP product label can also repeat when the same part occurs twice.
Use instance identity or independently verified placement to resolve
duplicate labels. For configured occurrences, preserve the part ID and full
configuration when available. An XML UID collision alone leaves the mapping
ambiguous. Zero XML mass and centroid values from unassigned materials
provide no reference. Both the nested and configurable examples had this
limitation.

A `403` on one candidate blocked the export with the current account.
Check document access and export permissions. Another large export was
interrupted on its first attempt and succeeded on a one-off retry with a
longer timeout.

Configured child assemblies, cross-document links, and the UI bulk STEP
comparison for the nested example remain open. Annual API
quota accounting for OAuth exports has not been verified.

[MathWorks exporter reference](https://www.mathworks.com/help/sm/ref/smexportonshape.html) ·
[MathWorks Onshape authentication and import guide](https://www.mathworks.com/help/sm/ug/onshape-import.html)
