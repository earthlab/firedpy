# Contributing to firedpy

Thank you for your interest in contributing to firedpy! This document provides guidelines for both internal Earth Lab team members and external contributors.

---

## Prerequisites
- Python 3.8 or higher
- Git

---

## Development Workflows

### For Internal Earth Lab / firedpy Team Members

#### Getting Started
1. **Contact the team**: Reach out to Nate, Aashish, or Adam to request a development branch
2. **Branch creation**: A dev branch will be created for you, cloned from `main`
3. **Clone your dev branch**: `git clone https://github.com/earthlab/firedpy.git -b your-dev-branch`
4. **Set up environment**: Follow the set-up steps in the [README.md](README.md)

#### Development Process
1. **Work on your `[your-initials]-dev` branch**: Make all changes on your assigned development branch
2. **Create pull requests**: PR your changes into `nonprod` branch
3. **Testing**: Ensure all tests pass before requesting review
4. **Review process**: Adam or Aashish will review your PRs

#### Branch Naming Convention
- Development branches: `[your-intials]-dev`
- Example: `NH-dev`

#### Branch Workflow Summary
```
[your-initials]-dev (your development branch)
  ↓
nonprod (staging/testing)
  ↓
main (production)
```

### For External Contributors (do we want to do this? better way to set up?)

#### Getting Started
1. **Fork the repository**: Click "Fork" on the GitHub repository page
2. **Clone your fork**: `git clone https://github.com/your-username/firedpy.git`
3. **Set up environment**: Follow the development setup steps below
4. **Create a feature branch**: `git checkout -b feature/your-feature-name`

#### Development Process
1. **Work on your feature branch**: Make changes on your feature branch
2. **Create pull requests**: PR your changes into the `nonprod` branch
3. **Review process**: Internal team members will review your PRs
4. **Testing**: Ensure all tests pass before requesting review

---

## How to Contribute

### Issues

#### **Issue Templates**
We use GitHub issue templates for reporting:

- **[Bug Report Template](.github/ISSUE_TEMPLATE/bug_report.md)** - For reporting bugs with scientific context
- **[Question Template](.github/ISSUE_TEMPLATE/question.md)** - For usage help and methodology questions

#### **Issue Review Process**
During each meeting, we'll review open issues and assign the following flags and attributes.

#### **Flags / Labels** (chat came up with these - do we like?)
| Label | Meaning |
|--------|----------|
| `prio:high` | Urgent or core algorithm changes (e.g., Fire Speed integration) |
| `prio:medium` | Secondary features or refactors |
| `prio:low` | Documentation, small cleanups |
| `type:bug` | Errors, failing tests, incorrect outputs |
| `type:feat` | New feature requests |
| `type:docs` | Documentation improvements |
| `type:perf` | Performance optimizations |
| `type:infra` | Build, branch organization, or repo maintenance |
| `area:io` | Data ingestion and preprocessing |
| `area:metrics` | Fire metrics (speed, FRP, etc.) |
| `area:agg` | Aggregation and spatial/temporal grouping |
| `area:docs` | User or API documentation |
| `good first issue` | Suitable for onboarding new contributors |

---

## Enhancements

#### **Enhancement Requests**
Use the **[Feature Request Template](.github/ISSUE_TEMPLATE/feature_request.md)** for suggesting new functionality.
See the [DEVELOPMENT.md](DEVELOPMENT.md) for a list of current development priorities.

---

## Pull Requests

#### **PR Templates**
We use different PR templates based on the scope of changes:

- **[Major Change/Release Template](.github/pull_request_template_major.md)** - For major features, breaking changes, or new releases
- **[Minor Fix/Update Template](.github/pull_request_template_minor.md)** - For small bug fixes and documentation updates

#### **PR Workflow**
- **Internal team**: PR your `[your-initials]-dev` branch into `nonprod`
- **External contributors**: PR your feature branch into `nonprod`
- Use template based on change scope
- Tag relevant issues in PR description
- Request review from maintainers (Adam or Aashish)

---

## Assigning Issues
- Each issue should have **one primary owner** and optional collaborators.
- Use GitHub **Projects** (Kanban view) to track issue progress:
  - `To Do`
  - `In Progress`
  - `Done`

---

## Coding Standards
- Follow PEP 8 style guidelines
- Use meaningful variable and function names
- Annotate your code!
- Write tests for new functionality

---

## Documentation
- Update docstrings for any new functions
- Add usage examples for new features
- Update README if needed

---

## Scientific Guidelines

### Data Validation
- Test changes + document
- Document any assumptions or limitations

### Reproducibility
- Use version control for all changes
- Document parameter choices
- Provide example workflows
- Consider backward compatibility

---

## Recognition

Contributors will be acknowledged in the project documentation and release notes. Significant contributions may warrant co-authorship on relevant publications.

Thank you for contributing to firedpy!
