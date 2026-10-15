// Independent structural schema check. DrawingML math is an Office extension;
// the SDK alone does not certify its nested payload or actual Office import.
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using DocumentFormat.OpenXml;
using DocumentFormat.OpenXml.Packaging;
using DocumentFormat.OpenXml.Validation;

if (args.Length != 2) throw new ArgumentException("usage: schema-tool corpus.json output-directory");
XNamespace math = "http://schemas.openxmlformats.org/officeDocument/2006/math";
XNamespace drawing = "http://schemas.openxmlformats.org/drawingml/2006/main";
const string a14 = "http://schemas.microsoft.com/office/drawing/2010/main";
var reportPath = Path.Combine(args[1], "result.json");
var report = JsonNode.Parse(File.ReadAllText(reportPath))!.AsObject();
var corpus = JsonNode.Parse(File.ReadAllText(args[0]))!.AsArray();
int fragments = 0;
int drawingPropertyLimitations = 0;
var validator = new OpenXmlValidator(FileFormatVersions.Office2010);
try
{
    foreach (var item in corpus)
    {
        if (item!["coverage"]?.GetValue<string>() == "excluded" || item["ooxml_rejected"] is not null) continue;
        var root = XElement.Load(Path.Combine(args[1], item["id"]!.GetValue<string>() + ".pptxmath"));
        if (root.Name != math + "oMathPara" || root.Descendants().Any(e => e.Name.Namespace != math && e.Name.Namespace != drawing))
            throw new InvalidDataException("Unexpected PowerPoint math namespace");
        foreach (var run in root.Descendants(math + "r"))
        {
            var names = run.Elements().Select(e => e.Name).ToArray();
            if (!names.SequenceEqual(new[] { math + "rPr", drawing + "rPr", math + "t" }) &&
                !names.SequenceEqual(new[] { drawing + "rPr", math + "t" }))
                throw new InvalidDataException("PowerPoint math run must use a:rPr/m:t");
        }
        foreach (var props in root.Descendants(drawing + "rPr"))
        {
            var size = (string?)props.Attribute("sz");
            var run = props.Parent;
            var argument = run?.Parent;
            var box = argument?.Parent;
            bool fixedToken = run?.Name == math + "r" && argument?.Name == math + "e" &&
                box?.Name == math + "box" && argument.Elements().Count() == 1 &&
                (string?)box.Element(math + "boxPr")?.Element(math + "ctrlPr")?.Element(drawing + "rPr")?.Attribute("sz") == "1800";
            if ((size != "1800" && !(fixedToken && new[] {"2160", "3240", "4320", "5400"}.Contains(size))) ||
                (string?)props.Element(drawing + "latin")?.Attribute("typeface") != "Cambria Math" ||
                !Regex.IsMatch((string?)props.Element(drawing + "solidFill")?.Element(drawing + "srgbClr")?.Attribute("val") ?? "", "^[0-9a-fA-F]{6}$"))
                throw new InvalidDataException("Unexpected native math font/size/colour profile");
        }
        var actualSizes = root.Descendants(drawing + "rPr").Select(p => (string?)p.Attribute("sz"))
            .Where(s => s != "1800").Select(s => int.Parse(s!)).ToArray();
        var expectedSizes = item["fixed_delimiter_sizes"]?.AsArray().Select(s => s!.GetValue<int>()).ToArray() ?? Array.Empty<int>();
        if (!actualSizes.SequenceEqual(expectedSizes))
            throw new InvalidDataException("Fixed delimiter sizes differ from the independent corpus profile");
        // Validate the shipped fragment without renaming text or properties.
        // The ISO SDK does not describe DrawingML rPr inside Office math;
        // those exact properties are checked against the strict profile above.
        var paragraph = new DocumentFormat.OpenXml.Math.Paragraph(root.ToString(SaveOptions.DisableFormatting));
        var errors = validator.Validate(paragraph).ToArray();
        var structuralErrors = errors.Where(e =>
            !(e.Description.Contains(drawing.NamespaceName + ":rPr") ||
              (e.Node?.NamespaceUri == drawing.NamespaceName && e.Node.LocalName == "rPr"))).ToArray();
        drawingPropertyLimitations += errors.Length - structuralErrors.Length;
        if (structuralErrors.Length != 0)
            throw new InvalidDataException(item["id"] + ": " + string.Join("; ", structuralErrors.Select(e => e.Description)));
        fragments++;
    }
    using var deck = PresentationDocument.Open(Path.Combine(args[1], "native-equations.pptx"), false);
    var packageErrors = validator.Validate(deck).ToArray();
    var outerErrors = packageErrors.Where(e =>
        !(e.Path?.XPath?.Contains("a14:m") ?? false) &&
        !(e.Node?.Ancestors().Any(n => n.NamespaceUri == a14 && n.LocalName == "m") ?? false)
    ).ToArray();
    if (outerErrors.Length != 0)
        throw new InvalidDataException(string.Join("; ", outerErrors.Select(e => e.Description)));
    report["schema_validation"] = new JsonObject
    {
        ["math_structure_fragments"] = fragments,
        ["drawing_run_property_sdk_limitations"] = drawingPropertyLimitations,
        ["outer_package_errors"] = outerErrors.Length,
        ["a14_sdk_limitation_errors"] = packageErrors.Length - outerErrors.Length,
        ["scope"] = "Unmodified shipped math fragments; declared DrawingML rPr SDK limitations; outer package; strict DrawingML run profile",
        ["full_a14_schema_certification"] = false,
        ["powerpoint_import_certified"] = false,
    };
}
catch (Exception error)
{
    report["passed"] = false;
    report["schema_validation"] = new JsonObject { ["error"] = error.Message };
    throw;
}
finally
{
    File.WriteAllText(reportPath, report.ToJsonString(new() { WriteIndented = true }));
}
Console.WriteLine($"Schema/profile checks passed for {fragments} math fragments; PowerPoint import remains a separate gate.");
