"""
Cubit C++ SDK customization knowledge base.

Covers building C++ plugins for Coreform Cubit that go beyond what the
in-process Python / PySide6 path can reach: new menu bar entries with
hotkeys, custom command panels with their own navigation-node tree, and
new command-line verbs that integrate with Cubit's command parser.

Sources:
  - Coreform webinar "Advanced customization in Coreform Cubit using the
    C++ SDK"  https://www.youtube.com/watch?v=u7GIqTKId98  (Carl McKelvey)
  - Reference plugin: Norbert Huber's "calculus" plugin (open source,
    tied directly into Cubit -- a complete worked example).
  - Coreform forum: https://forum.coreform.com  (Norbert is active there).

Scope split with sibling knowledge modules:
  - `custom_toolbar.py` -- in-process Python / PySide6 toolbars + dialogs.
    First choice for most customization; ~no build system, no Qt6 install.
  - this module        -- C++ SDK plugins.  Only when Python cannot do it:
    + Hotkey-bound NEW custom commands.
    + Custom command panels integrated into Cubit's navigation-node tree
      (the left-side command panel area).
    + Deep menu/toolbar manipulation beyond Tools->Custom Toolbar Editor.
    + Sharing a complex C++ workflow as one binary across an org.
  - `panel_conventions.py` -- Radia-NGSolve external Python 3.12 panels.
    Unrelated to this module; lives outside Cubit's process.
"""

CPP_SDK_OVERVIEW = """
# Cubit C++ SDK Plugin Overview

The C++ SDK lets you build native plugins (`.dll` on Windows, `.so` on
Linux) that load into Cubit at startup and extend the GUI / command set.
Unlike the in-process Python toolbar workflow, the SDK exposes Claro
(Cubit's GUI wrapper) directly and gives you the same hooks the
built-in Cubit modules use.

## When the C++ SDK is the right choice

Reach for the SDK when you need any of these:

- **Hotkey-bound new commands**.  Qt supports `&M` syntax (Alt+M) on
  menu actions you ADD, but existing commands cannot have hotkeys
  retrofitted from Python.  A new C++ command panel + action can.
- **Custom command panels** integrated into the left-side
  navigation-node tree (Geometry / Mesh / Boundary Conditions / etc.).
  Python toolbars can only spawn floating dialogs.
- **New command-line verbs** that integrate with Cubit's command
  parser (so users can `play` journal files that drive them).
- **Org-wide distributable binary** for proprietary workflows where
  Python source would be too exposed (Python files in Cubit toolbars
  are always plain text; obfuscation is "more work than it's worth").

## When the Python toolbar path is enough -- and preferred

- One-off parameter dialogs.
- Multi-step workflow encapsulation (clean -> imprint -> mesh -> export
  as one button).
- Anything an LLM can generate as a single PySide6 file.

See `custom_toolbar.py` for the Python path.  Don't reach for C++ SDK
when a 50-line Python script does the job.

## Reference open-source plugin

The single best worked example is **Norbert Huber's "calculus" plugin**:
a full C++ SDK plugin tied directly into Cubit, available on the
Coreform forum.  Read its CMakeLists.txt + menu_manager + panel
factory before writing your own.

## License / packaging notes

- The Cubit SDK headers ship inside the Cubit install (`<install>/SDK/`)
  and are intended to be COPIED into your project tree; do not modify
  in place.
- Qt 6.8 is the upstream version Cubit links against.  Use the LGPL Qt
  distribution unless you have a commercial Qt license; LGPL imposes
  the usual dynamic-linking constraints on your plugin.
- Plugins are loaded at Cubit startup; users must **restart Cubit**
  after dropping a new plugin folder into Tools -> Plugins -> Add.
  There is no hot-reload.
"""

CPP_SDK_PREREQUISITES = """
# Cubit C++ SDK Prerequisites

| Tool        | Version             | Reason                              |
|-------------|---------------------|-------------------------------------|
| Cubit       | 2025.12+ recommended| SDK directory shipped + supported   |
| CMake       | 3.20+               | Cross-platform build system Cubit uses |
| C++ compiler| C++17 minimum       | Cubit headers require it            |
|   Windows   | Visual Studio 2022  | MSVC 19.30+ (vcvarsx64 environment) |
|   Linux     | GCC 13+ or Clang 16+| Either works                        |
| Qt          | 6.8 (LGPL OK)       | Cubit links against 6.8 exactly     |

## Locating the SDK

The SDK directory ships inside the Cubit install:

```
Windows:  C:\\Program Files\\Coreform Cubit 2025.12\\SDK\\
Linux:    /opt/Coreform-Cubit-2025.12/SDK/
```

Copy this entire `SDK/` tree to your own development area so you can
modify CMakeLists.txt without touching the install.  The Cubit install
itself is typically read-only.

## CMakeLists.txt hints

The SDK's top-level `CMakeLists.txt` needs to find Cubit and Qt6.
The canonical pattern (taken from the components example) is:

```cmake
set(Cubit_DIR "C:/Program Files/Coreform Cubit 2025.12/cmake"
    CACHE PATH "Cubit cmake config dir")
set(Qt6_DIR   "C:/Qt/6.8.0/msvc2022_64/lib/cmake/Qt6"
    CACHE PATH "Qt6 cmake config dir")

find_package(Cubit REQUIRED)
find_package(Qt6 6.8 REQUIRED COMPONENTS Core Widgets)
```

Setting these as `CACHE PATH` (not plain `set()`) means a developer can
override them from the cmake command line without editing the file.

## Build directory layout

Cubit's convention is to build OUT-OF-SOURCE, parallel to the source
tree (not inside it):

```
my-plugin/
  src/
  CMakeLists.txt          <-- source (version-controlled)
build-my-plugin/          <-- parallel build dir (gitignored)
  Debug/   my_plugin.dll
  Release/ my_plugin.dll
```

From `build-my-plugin/`:

```
# Windows (from an x64 Native Tools Command Prompt):
cmake -G "Visual Studio 17 2022" -A x64 ../my-plugin
cmake --build . --config Release

# Linux:
cmake -G Ninja ../my-plugin
cmake --build .
```

## Mixing Debug/Release on Windows

**Do not mix MSVC Debug and Release builds.**  If you build the plugin
as Debug but Cubit ships as Release (it does), loading the plugin will
either fail silently or crash the GUI on first dialog open.  Always
build your plugin as **Release** unless you also have a Debug-built
Cubit (rare).
"""

CPP_SDK_INSTALL_AND_LOAD = """
# Loading a Cubit C++ SDK Plugin

After `cmake --build . --config Release` produces
`build-my-plugin/Release/my_plugin.dll`:

## Step 1: Load via Tools -> Plugins

In the Cubit GUI:

1. **Tools menu -> Plugins**
2. Click **Add**
3. Select the **folder** containing the .dll (NOT the .dll itself --
   Cubit scans the folder for plugin binaries).
   E.g.: `C:/dev/build-my-plugin/Release/`
4. Close the Plugins dialog.

## Step 2: Restart Cubit

Plugins are loaded **only at Cubit startup**.  After Add, you MUST exit
and relaunch Cubit for the plugin to appear.  This is the single most
common cause of "I added the plugin but nothing happened" -- the user
forgot to restart.

## Step 3: Verify

After restart, reopen Tools -> Plugins; your plugin name should be in
the list.  If your plugin adds a menu / toolbar / panel, those should
also be visible now.

## Distribution to other machines

Two paths:

1. **Per-user**: each user copies the plugin folder somewhere on their
   machine and does Tools -> Plugins -> Add themselves.
2. **Org-wide**: install the plugin folder into a shared / network
   location (or a global Cubit plugin dir) and tell users to add it.
   Cubit does not have a system-wide "auto-load" mechanism; each user
   still does the Add+restart once.

Note that Cubit's "Plugins" dialog remembers the folder path, so a
user only needs to Add once per machine; subsequent Cubit restarts
auto-load the plugin from that path.

## Common failure modes

| Symptom                              | Likely cause                       |
|--------------------------------------|------------------------------------|
| Plugin not in Tools->Plugins after restart | Wrong folder selected, or .dll not actually built (check filesystem) |
| Cubit crashes on startup            | Debug/Release ABI mismatch -- rebuild plugin in Release |
| Plugin loads but menu/panel missing | Missing `Q_OBJECT` macro on a class, or moc/uic step skipped in CMake |
| Hotkey doesn't trigger              | Forgot the `&` in the menu action text (`"&MyAction"` for Alt+M) |
"""

CPP_SDK_COMPONENTS_EXAMPLE = """
# The SDK `components` Example -- Anatomy

Cubit ships with an example plugin under `SDK/components/`.  This is the
single best starting point; it demonstrates ALL the customization hooks
in one buildable project.

## File layout

```
SDK/components/
  CMakeLists.txt              -- build (add Cubit/Qt6 hints here)
  component.cpp / .h          -- plugin entry point, lifecycle hooks
  menu_manager.cpp / .h       -- add menu bar items + hotkeys
  toolbar_manager.cpp / .h    -- add toolbar buttons
  panel_factory.cpp / .h      -- create custom command panels
  navigation_nodes.xml        -- declare panels' location in the left-side tree
  panels/
    my_panel.cpp / .h         -- one C++ class per panel
    my_panel.ui               -- Qt Designer .ui file (optional)
```

## Plugin entry point

The `component.cpp` file exports the plugin entry hooks Cubit calls at
load time.  It typically:

1. Asks Claro (Cubit's GUI wrapper) for the main-window instance.
2. Calls `menu_manager::setup_custom_menu(claro)`.
3. Calls `toolbar_manager::setup_custom_toolbar(claro)`.
4. Calls `panel_factory::register_panels(claro)`.

There's also a `s_custom_initialize` static flag pattern used to
guard against re-initialization if Cubit calls the entry hook twice:

```cpp
static bool s_custom_initialize = false;
void component_init() {
    if (s_custom_initialize) return;
    auto claro = get_claro_instance();
    if (!claro) return;
    // ... setup menus, toolbars, panels ...
    s_custom_initialize = true;
}
```

## Why the design splits Claro from Cubit core

Carl McKelvey's recap (~20 years ago):

> "We decided we want the user interface and the core code to be
>  separate.  So we added in interfaces inside of Cubit that are
>  exposed -- we can pass commands down through them, we can query
>  things back from Cubit.  So anything that you can see or do in the
>  Cubit GUI actually goes through an interface to get into Cubit."

Implication: the Python API and the C++ SDK both ride the same Claro
interface.  The C++ SDK simply gets more of it (private Claro
methods, panel registration, command parser hooks) than the Python
wrapper exposes.

## Customizing the menu bar

`menu_manager.cpp` shows the canonical pattern:

```cpp
void setup_custom_menu(ClaroInterface* claro) {
    std::vector<QAction*> actions;

    auto* something_cool = new QAction(QObject::tr("&MyAction"), nullptr);
    QObject::connect(something_cool, &QAction::triggered, []() {
        // do something cool
        claro->command_window()->print("Custom command fired\\n");
    });
    actions.push_back(something_cool);

    claro->add_menu("My Component", actions);
}
```

Key points:
- `&` in the action text marks the next character as the Alt-key
  hotkey (`&M` -> Alt+M).
- `QAction::triggered` is connected to a slot (lambda OK).
- Inside the slot, use `claro->command_window()->print(...)` to log to
  Cubit's bottom command window.
- The action shows up in a new top-level menu named "My Component".

## Customizing toolbars

`toolbar_manager.cpp` follows the same pattern:

```cpp
auto* tb = claro->add_toolbar("MyToolbar");
auto* btn = new QAction(QIcon(":/icons/my.png"), QObject::tr("Top View"), nullptr);
QObject::connect(btn, &QAction::triggered, []() { /* ... */ });
tb->addAction(btn);
```

Icons can come from a Qt resource file (`.qrc`) compiled into the
plugin via `qt6_add_resources` in CMakeLists.txt.

## Custom command panels (the hard part)

Carl's summary in the webinar:

> "Setting up a custom command panel is admittedly more complex."

The four pieces you have to author:

1. **The Qt widget itself** (the panel UI).  Code it directly or use
   the Qt Designer `.ui` file.  Most of the built-in Cubit panels were
   designed in Qt Designer.
2. **Navigation nodes** -- the tree on the left side of the Cubit GUI
   (Geometry / Mesh / etc.).  Your panels need to be attached to a
   node.  Two ways to declare:
     - Programmatically, from C++ (see the components example).
     - Via an XML file (also in the components example).  The XML path
       requires installing the XML where Cubit runs from, which means
       root/admin privilege on the install location -- programmatic
       registration is usually less painful.
3. **A factory** that instantiates the panel widget on demand and
   associates it with the navigation node.  This is the
   `panel_factory.cpp` pattern.
4. **Wiring back to Cubit commands** -- the panel's "Apply" button
   typically builds a Cubit command string and submits it via
   `claro->command_interface()->run_command(...)`.

If you only need a floating dialog (not docked into the navigation
tree), the Python toolbar path is dramatically simpler.  Reach for
custom command panels only when the navigation-tree integration is
required.
"""

CPP_SDK_CUSTOM_COMMANDS = """
# Adding New Cubit Command-Line Verbs (C++ SDK)

The C++ SDK lets you register entirely new command verbs that show up
in Cubit's command parser.  Once registered, users can type your
command at the Cubit command prompt, or `play` a journal file that
invokes it -- exactly like a built-in command.

This is a capability the in-process Python toolbar workflow cannot
match: Python toolbars only run callbacks; they don't extend the
command grammar.

## When to add a new command

- The action should be `play`-replayable from a `.jou` file.
- The action should be hotkey-bindable (only NEW commands can get
  hotkeys; existing commands cannot have hotkeys retrofitted -- see
  webinar QA).
- The action benefits from being scripted in a journal file rather
  than triggered from a GUI button.

## Registration sketch

The SDK example shows the `command_factory` pattern.  Conceptually:

```cpp
class MyCommand : public CubitCommand {
public:
    virtual CubitStatus execute() override {
        // parse parameters from the command line
        // run the action
        // return CUBIT_SUCCESS or CUBIT_FAILURE
    }
    virtual std::string syntax() const override { return "my_command <int>"; }
    virtual std::string help() const override { return "..."; }
};

void register_commands(ClaroInterface* claro) {
    claro->command_parser()->register_command("my_command", []() {
        return std::make_unique<MyCommand>();
    });
}
```

After plugin load, users can do:

```
cubit> my_command 42
cubit> play my_journal.jou       # if the .jou contains `my_command 42`
```

## Hotkey for the new command

```cpp
auto* action = new QAction(QObject::tr("&MyCommand"), nullptr);
action->setShortcut(QKeySequence("Ctrl+Shift+M"));
QObject::connect(action, &QAction::triggered, []() {
    claro->command_interface()->run_command("my_command 42");
});
```

The `&M` Alt-hotkey works on menu items; `setShortcut` with
`QKeySequence("Ctrl+Shift+M")` works for arbitrary key combinations.
"""

CPP_SDK_TROUBLESHOOTING = """
# Cubit C++ SDK Troubleshooting

## "Plugin doesn't show up after Tools->Plugins->Add"

- Did you select the FOLDER containing the .dll, not the .dll itself?
- Did you exit and relaunch Cubit?  Plugins load only at startup.
- Check `cmake --build` actually produced a .dll in the folder you
  added.  An empty Release/ folder gives no feedback.

## "Cubit crashes on startup after adding my plugin"

- Almost always Debug/Release ABI mismatch on Windows.  Rebuild your
  plugin in Release.  Don't link Debug Qt or Debug MSVCRT against a
  Release Cubit.
- Second-most-common: Qt version drift.  Cubit links Qt 6.8 exactly.
  Building your plugin against Qt 6.5 or 6.9 can crash on first widget
  creation if any class has shifted vtable layout.

## "MOC / UIC didn't run -- missing virtual methods at link time"

CMake's auto-MOC needs to be turned on explicitly:

```cmake
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
```

Without these, MOC is skipped silently and you get link errors for
`Q_OBJECT`-decorated classes.

## "Menu item appears but hotkey doesn't trigger"

- Forgot the `&` in the action text.  `"My Action"` -> no hotkey;
  `"&My Action"` -> Alt+M.
- The hotkey letter must be unique among visible menu items at the
  same level; Qt will silently disable a duplicate.

## "Panel appears but is empty / has no buttons"

- Check that the `.ui` file (if you use Qt Designer) is being run
  through UIC at build time.  `set(CMAKE_AUTOUIC ON)` plus
  `target_sources(my_plugin PRIVATE my_panel.ui)` is the standard fix.
- Check that the panel's `Q_OBJECT` is present and the class is
  exported (not just defined locally inside the .cpp).

## "Can I see what's happening inside Cubit at load time?"

Cubit's command window (bottom of the GUI) shows plugin load
diagnostics.  For deeper inspection, attach Visual Studio's debugger
to the Cubit process after launch and load symbols for your .dll.

## Where to get help

Coreform forum: https://forum.coreform.com/

Norbert Huber (calculus plugin author) is active there and answers
SDK-specific questions.  Post a minimal reproduction (your
CMakeLists.txt + the smallest .cpp that fails) for fastest response.
"""


CPP_SDK_TOPICS = {
    "overview": CPP_SDK_OVERVIEW,
    "prerequisites": CPP_SDK_PREREQUISITES,
    "install_and_load": CPP_SDK_INSTALL_AND_LOAD,
    "components_example": CPP_SDK_COMPONENTS_EXAMPLE,
    "custom_commands": CPP_SDK_CUSTOM_COMMANDS,
    "troubleshooting": CPP_SDK_TROUBLESHOOTING,
}

# Aliases for natural-language access from the MCP server.
CPP_SDK_ALIASES = {
    "intro": "overview",
    "when_to_use": "overview",
    "setup": "prerequisites",
    "cmake": "prerequisites",
    "qt": "prerequisites",
    "qt6": "prerequisites",
    "build": "prerequisites",
    "compile": "prerequisites",
    "load": "install_and_load",
    "install": "install_and_load",
    "plugins_dialog": "install_and_load",
    "restart": "install_and_load",
    "components": "components_example",
    "example": "components_example",
    "claro": "components_example",
    "menu": "components_example",
    "menu_manager": "components_example",
    "toolbar": "components_example",
    "panel": "components_example",
    "panel_factory": "components_example",
    "navigation_nodes": "components_example",
    "commands": "custom_commands",
    "command_verb": "custom_commands",
    "hotkey": "custom_commands",
    "shortcut": "custom_commands",
    "debug": "troubleshooting",
    "crash": "troubleshooting",
    "moc": "troubleshooting",
    "uic": "troubleshooting",
}


def get_cpp_sdk_documentation(topic: str = "overview") -> str:
    """Return the C++ SDK knowledge body for a topic.

    Args:
        topic: One of CPP_SDK_TOPICS keys, an alias in CPP_SDK_ALIASES,
            or "all" to get every section concatenated.

    Returns:
        Markdown documentation string, or an error message naming
        the available topics if the input is unknown.
    """
    topic = (topic or "overview").lower().strip()
    if topic == "all":
        order = ["overview", "prerequisites", "install_and_load",
                 "components_example", "custom_commands",
                 "troubleshooting"]
        return "\n\n".join(CPP_SDK_TOPICS[k] for k in order)
    resolved = CPP_SDK_ALIASES.get(topic, topic)
    if resolved in CPP_SDK_TOPICS:
        return CPP_SDK_TOPICS[resolved]
    valid = sorted(CPP_SDK_TOPICS.keys())
    return (
        f"Unknown cpp_sdk topic: {topic!r}.  "
        f"Valid topics: {valid}.  "
        f"Aliases: {sorted(CPP_SDK_ALIASES.keys())}."
    )
