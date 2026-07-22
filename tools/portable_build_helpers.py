from pathlib import Path


def discover_python_modules(site_packages: Path, package_name: str) -> list[str]:
    """List import names without relying on namespace-package recursion."""
    package_root = Path(site_packages) / package_name
    modules: set[str] = set()
    for source in package_root.rglob("*.py"):
        relative = source.relative_to(package_root)
        parts = list(relative.with_suffix("").parts)
        if parts[-1] == "__init__":
            parts.pop()
        modules.add(".".join([package_name, *parts]).rstrip("."))
    return sorted(modules)
