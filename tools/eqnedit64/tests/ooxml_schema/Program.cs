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
XNamespace word = "http://schemas.openxmlformats.org/wordprocessingml/2006/main";
const string a14 = "http://schemas.microsoft.com/office/drawing/2010/main";
var reportPath = Path.Combine(args[1], "result.json");
var report = JsonNode.Parse(File.ReadAllText(reportPath))!.AsObject();
var corpus = JsonNode.Parse(File.ReadAllText(args[0]))!.AsArray();
int fragments = 0;
var validator = new OpenXmlValidator(FileFormatVersions.Office2010);
try
{
    foreach (var item in corpus)
    {
        if (item!["coverage"]?.GetValue<string>() == "excluded") continue;
        var root = XElement.Load(Path.Combine(args[1], item["id"]!.GetValue<string>() + ".pptxmath"));
        if (root.Name != math + "oMathPara" || root.Descendants().Any(e => e.Name.Namespace != math && e.Name.Namespace != drawing))
            throw new InvalidDataException("Unexpected PowerPoint math namespace");
        foreach (var run in root.Descendants(math + "r"))
        {
            var names = run.Elements().Select(e => e.Name).ToArray();
            if (!names.SequenceEqual(new[] { math + "rPr", drawing + "rPr", drawing + "t" }) &&
                !names.SequenceEqual(new[] { drawing + "rPr", drawing + "t" }))
                throw new InvalidDataException("PowerPoint math run must use a:rPr/a:t");
        }
        foreach (var props in root.Descendants(drawing + "rPr"))
        {
            if ((string?)props.Attribute("sz") != "1800" ||
                (string?)props.Element(drawing + "latin")?.Attribute("typeface") != "Cambria Math" ||
                !Regex.IsMatch((string?)props.Element(drawing + "solidFill")?.Element(drawing + "srgbClr")?.Attribute("val") ?? "", "^[0-9a-fA-F]{6}$"))
                throw new InvalidDataException("Unexpected native math font/size/colour profile");
        }
        // Validation-only copy: remove DrawingML run formatting and normalize
        // text tags so the independent ISO math-structure schema can check
        // fractions/scripts/nary/matrices/functions. This is not a production
        // translator and is not claimed to validate the a14 extension itself.
        var structural = new XElement(root);
        foreach (var props in structural.Descendants(drawing + "rPr").ToList())
            props.ReplaceWith(new XElement(word + "rPr"));
        foreach (var text in structural.Descendants(drawing + "t")) text.Name = math + "t";
        var paragraph = new DocumentFormat.OpenXml.Math.Paragraph(structural.ToString(SaveOptions.DisableFormatting));
        var errors = validator.Validate(paragraph).ToArray();
        if (errors.Length != 0)
            throw new InvalidDataException(item["id"] + ": " + string.Join("; ", errors.Select(e => e.Description)));
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
        ["outer_package_errors"] = outerErrors.Length,
        ["a14_sdk_limitation_errors"] = packageErrors.Length - outerErrors.Length,
        ["scope"] = "ISO math structure after validation-only run normalization; outer package; strict DrawingML run profile",
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
